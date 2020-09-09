# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.analyzer.defect_energy import make_defect_energies
from pydefect.analyzer.defect_energy_plotter import DefectEnergyMplPlotter
from pydefect.analyzer.defect_structure_analyzer import DefectStructureAnalyzer
from pydefect.analyzer.eigenvalue_plotter import EigenvalueMplPlotter
from pydefect.analyzer.make_band_edge_state import make_band_edge_state
from pydefect.analyzer.make_defect_energy import make_single_defect_energy
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo, \
    CompositionEnergy, replace_comp_energy
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiagMpl2DMplPlotter, \
    ChemPotDiagMpl3DMplPlotter
from pydefect.chem_pot_diag.make_chem_pot_diag import make_chem_pot_diag_from_mp
from pydefect.cli.main_tools import sanitize_matrix
from pydefect.cli.vasp.make_band_edge_eigenvalues import \
    make_band_edge_eigenvalues
from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp
from pydefect.cli.vasp.make_edge_characters import MakeEdgeCharacters
from pydefect.cli.vasp.make_efnv_correction import \
    make_efnv_correction
from pydefect.cli.vasp.make_poscars_from_query import make_poscars_from_query
from pydefect.cli.vasp.make_unitcell import make_unitcell_from_vasp
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialMplPlotter
from pydefect.defaults import defaults
from pydefect.input_maker.add_interstitial import append_interstitial
from pydefect.input_maker.defect_entries_maker import DefectEntriesMaker
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_set_maker import DefectSetMaker
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.input_maker.supercell_maker import SupercellMaker
from pydefect.util.error_classes import CpdNotSupportedError
from pydefect.util.mp_tools import MpQuery
from pymatgen.io.vasp import Vasprun, Outcar, Procar
from pymatgen.util.string import latexify
from vise.util.logger import get_logger

logger = get_logger(__name__)


def print_file(args):
    print(args.obj)


def make_unitcell(args):
    unitcell = make_unitcell_from_vasp(vasprun_band=args.vasprun_band,
                                       outcar_band=args.outcar_band,
                                       outcar_dielectric=args.outcar_dielectric)
    unitcell.to_json_file()


def make_competing_phase_dirs(args):
    query = MpQuery(element_list=args.elements, e_above_hull=args.e_above_hull)
    make_poscars_from_query(materials_query=query.materials, path=Path.cwd())


def make_chem_pot_diag(args) -> None:
    if args.elements:
        cpd = make_chem_pot_diag_from_mp(additional_elements=args.elements,
                                         target=args.target,
                                         atom_energy_yaml=args.functional)
    else:
        comp_es = []
        for d in args.dirs:
            vasprun = Vasprun(d / defaults.vasprun)
            composition = vasprun.final_structure.composition
            energy = float(vasprun.final_energy)  # type is FloatWithUnit
            comp_es.append(CompositionEnergy(composition, energy, "local"))
        if args.update:
            cpd = ChemPotDiag.from_yaml(args.yaml)
            replace_comp_energy(cpd, comp_es)
        else:
            cpd = ChemPotDiag(comp_es, args.target)
    cpd.to_yaml(args.yaml)


def plot_chem_pot_diag(args) -> None:
    cpd = ChemPotDiag.from_yaml(args.yaml)
    if cpd.dim == 2:
        plotter = ChemPotDiagMpl2DMplPlotter(CpdPlotInfo(cpd))
    elif cpd.dim == 3:
        plotter = ChemPotDiagMpl3DMplPlotter(CpdPlotInfo(cpd))
    else:
        raise CpdNotSupportedError("Number of elements must be 2 or 3. "
                                   f"Now {cpd.vertex_elements}.")
    plt = plotter.draw_diagram()
    plt.savefig(fname="cpd.pdf")
    plt.show()


def make_supercell(args):
    if args.matrix:
        matrix = sanitize_matrix(args.matrix)
        maker = SupercellMaker(args.unitcell, matrix)
    else:
        kwargs = {}
        if args.min_num_atoms:
            kwargs["min_num_atoms"] = args.min_num_atoms
        if args.max_num_atoms:
            kwargs["max_num_atoms"] = args.max_num_atoms
        maker = SupercellMaker(args.unitcell, **kwargs)

    maker.supercell.structure.to(filename="SPOSCAR")
    maker.supercell_info.to_json_file()


def append_interstitial_to_supercell_info(args):
    supercell_info = append_interstitial(args.supercell_info,
                                         args.base_structure,
                                         args.frac_coords)
    supercell_info.to_json_file()


def pop_interstitial_from_supercell_info(args):
    supercell_info = args.supercell_info
    supercell_info.interstitials.pop(args.index - 1)
    supercell_info.to_json_file()


def make_defect_set(args):
    supercell_info = loadfn("supercell_info.json")
    if args.oxi_states:
        oxi_states = dict(zip(args.oxi_states[::2], args.oxi_states[1::2]))
    else:
        oxi_states = None
    maker = DefectSetMaker(supercell_info,
                           oxi_states,
                           args.dopants,
                           keywords=args.kwargs)
    maker.defect_set.to_yaml()


def make_defect_entries(args):
    supercell_info: SupercellInfo = loadfn("supercell_info.json")
    perfect = Path("perfect")

    try:
        perfect.mkdir()
        logger.info("Making perfect dir...")
        supercell_info.structure.to(filename=perfect / "POSCAR")
    except FileExistsError:
        logger.info(f"perfect dir exists, so skipped...")

    defect_set = DefectSet.from_yaml()
    maker = DefectEntriesMaker(supercell_info, defect_set)

    for defect_entry in maker.defect_entries:
        dir_path = Path(defect_entry.full_name)
        try:
            dir_path.mkdir()
            logger.info(f"Making {dir_path} dir...")
            defect_entry.perturbed_structure.to(filename=dir_path / "POSCAR")
            defect_entry.to_json_file(filename=dir_path / "defect_entry.json")
            defect_entry.to_prior_info(filename=dir_path / "prior_info.yaml")
        except FileExistsError:
            logger.info(f"{dir_path} dir exists, so skipped...")


def make_calc_results(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        calc_results = make_calc_results_from_vasp(
            vasprun=Vasprun(d / defaults.vasprun),
            outcar=Outcar(d / defaults.outcar))
        calc_results.to_json_file(filename=Path(d) / "calc_results.json")


def make_efnv_correction_from_vasp(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        defect_entry: DefectEntry = loadfn(d / "defect_entry.json")
        calc_results = loadfn(d / "calc_results.json")
        efnv = make_efnv_correction(defect_entry.charge,
                                    calc_results,
                                    args.perfect_calc_results,
                                    args.unitcell.dielectric_constant)
        efnv.to_json_file(d / "correction.json")

        title = defect_entry.full_name
        plotter = SitePotentialMplPlotter(title=title, efnv_correction=efnv)
        plotter.construct_plot()
        plotter.plt.savefig(fname=d / "correction.pdf")
        plotter.plt.clf()


def make_defect_eigenvalues(args):
    supercell_vbm = args.perfect_calc_results.vbm
    supercell_cbm = args.perfect_calc_results.cbm
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        defect_entry = loadfn(d / "defect_entry.json")
        title = defect_entry.name
        vasprun = Vasprun(d / defaults.vasprun)
        band_edge_eigenvalues = make_band_edge_eigenvalues(
            vasprun, supercell_vbm, supercell_cbm)
        band_edge_eigenvalues.to_json_file(d / "band_edge_eigenvalues.json")
        plotter = EigenvalueMplPlotter(
            title=title, band_edge_eigenvalues=band_edge_eigenvalues,
            supercell_vbm=supercell_vbm, supercell_cbm=supercell_cbm)
        plotter.construct_plot()
        plotter.plt.savefig(fname=d / "eigenvalues.pdf")
        plotter.plt.clf()


def make_edge_characters(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        vasprun = Vasprun(d / defaults.vasprun)
        procar = Procar(d / defaults.procar)
        outcar = Outcar(d / defaults.outcar)
        calc_results = loadfn(d / "calc_results.json")
        structure_analyzer = DefectStructureAnalyzer(
            calc_results.structure, args.perfect_calc_results.structure)
        edge_characters = MakeEdgeCharacters(
            procar, vasprun, outcar,
            structure_analyzer.neighboring_atom_indices).edge_characters
        edge_characters.to_json_file(d / "edge_characters.json")


def make_edge_states(args):
    for d in args.dirs:
        print(f"-- {d}")
        edge_states = []
        edge_characters = loadfn(d / "edge_characters.json")
        for spin, edge_character, ref in zip(["spin up  ", "spin down"],
                                             edge_characters,
                                             args.perfect_edge_characters):
            edge_state = make_band_edge_state(edge_character, ref)
            edge_states.append(edge_state)
            print(spin, edge_state)

        BandEdgeStates(edge_states).to_json_file(d / "band_edge_states.json")


def make_defect_formation_energy(args):
    title = latexify(
        args.perfect_calc_results.structure.composition.reduced_formula)
    chem_pot_diag = ChemPotDiag.from_yaml(args.chem_pot_diag)
    abs_chem_pot = chem_pot_diag.abs_chem_pot_dict(args.label)

    single_energies = []
    for d in args.dirs:
        if args.skip_shallow and loadfn(d / "band_edge_states.json").is_shallow:
            continue
        single_energies.append(
            make_single_defect_energy(args.perfect_calc_results,
                                      loadfn(d / "calc_results.json"),
                                      loadfn(d / "defect_entry.json"),
                                      abs_chem_pot,
                                      loadfn(d / "correction.json")))

    defect_energies = make_defect_energies(single_energies)
    if args.print:
        print("         charge          E_f   correction    ")
        for e in defect_energies:
            print(e)
            print("")

        print("-- cross points -- ")
        for e in defect_energies:
            print(e.name)
            print(e.cross_points(args.unitcell.vbm, args.unitcell.cbm))
            print("")
        return

    plotter = DefectEnergyMplPlotter(title=title,
                                     defect_energies=defect_energies,
                                     vbm=args.unitcell.vbm,
                                     cbm=args.unitcell.cbm,
                                     supercell_vbm=args.perfect_calc_results.vbm,
                                     supercell_cbm=args.perfect_calc_results.cbm,
                                     y_range=args.y_range)

    plotter.construct_plot()
    plotter.plt.savefig(f"energy_{args.label}.pdf")
