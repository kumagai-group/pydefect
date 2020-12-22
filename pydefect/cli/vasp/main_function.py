# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.analyzer.dash_components.cpd_energy_dash import \
    CpdEnergyComponent, make_energies
from pydefect.analyzer.defect_energy import slide_energy
from pydefect.analyzer.defect_energy_plotter import DefectEnergyMplPlotter
from pydefect.analyzer.defect_structure_analyzer import DefectStructureAnalyzer
from pydefect.analyzer.eigenvalue_plotter import EigenvalueMplPlotter
from pydefect.analyzer.make_band_edge_state import make_band_edge_state
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
from pydefect.cli.vasp.make_gkfo_correction import make_gkfo_correction
from pydefect.cli.vasp.make_poscars_from_query import make_poscars_from_query
from pydefect.cli.vasp.make_unitcell import make_unitcell_from_vasp
from pydefect.corrections.site_potential_plotter import \
    SitePotentialMplPlotter
from pydefect.defaults import defaults
from pydefect.input_maker.add_interstitial import append_interstitial
from pydefect.input_maker.defect_entries_maker import DefectEntriesMaker
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_set_maker import DefectSetMaker
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.input_maker.supercell_maker import SupercellMaker
from pydefect.util.mp_tools import MpQuery
from pymatgen.io.vasp import Vasprun, Outcar, Procar
from pymatgen.util.string import latexify
from vise.util.logger import get_logger

logger = get_logger(__name__)


def print_file(args):
    print(args.obj)


def make_unitcell(args):
    unitcell = make_unitcell_from_vasp(
        vasprun_band=args.vasprun_band,
        outcar_band=args.outcar_band,
        outcar_dielectric_clamped=args.outcar_dielectric_clamped,
        outcar_dielectric_ionic=args.outcar_dielectric_ionic)
    unitcell.to_json_file()


def make_competing_phase_dirs(args):
    query = MpQuery(element_list=args.elements, e_above_hull=args.e_above_hull)
    make_poscars_from_query(materials_query=query.materials, path=Path.cwd())


def make_chem_pot_diag(args) -> None:
    if args.elements:
        cpd = make_chem_pot_diag_from_mp(additional_elements=args.elements,
                                         target=args.target,
                                         atom_energy_yaml=args.atom_energy_yaml)
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
    if cpd.dim == 1:
        logger.warning("Single element is not supported for the plot.")
        return

    print(cpd)
    if cpd.dim == 2:
        plotter = ChemPotDiagMpl2DMplPlotter(CpdPlotInfo(cpd))
    elif cpd.dim == 3:
        plotter = ChemPotDiagMpl3DMplPlotter(CpdPlotInfo(cpd))
    else:
        logger.info("Number of elements must be 2 or 3. "
                    f"Now, elements are {cpd.vertex_elements}.")
        return
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
    supercell_info: SupercellInfo = append_interstitial(args.supercell_info,
                                                        args.base_structure,
                                                        args.frac_coords)
    supercell_info.to_json_file()


def pop_interstitial_from_supercell_info(args):
    supercell_info = args.supercell_info
    assert args.index > 0
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
        plotter = SitePotentialMplPlotter.from_efnv_corr(
            title=title, efnv_correction=efnv)
        plotter.construct_plot()
        plotter.plt.savefig(fname=d / "correction.pdf")
        plotter.plt.clf()


def make_gkfo_correction_from_vasp(args):
    gkfo = make_gkfo_correction(
        efnv_correction=args.initial_efnv_correction,
        additional_charge=args.charge_diff,
        final_calc_results=args.final_calc_results,
        initial_calc_results=args.initial_calc_results,
        diele_tensor=args.unitcell.dielectric_constant,
        ion_clamped_diele_tensor=args.unitcell.ele_dielectric_const)
    print(gkfo)
    gkfo.to_json_file("gkfo_correction.json")
    plotter = SitePotentialMplPlotter.from_gkfo_corr(
        title="GKFO correction", gkfo_correction=gkfo)
    plotter.construct_plot()
    plotter.plt.savefig(fname="gkfo_correction.pdf")
    plotter.plt.clf()


def make_defect_eigenvalues(args):
    supercell_vbm = args.perfect_calc_results.vbm
    supercell_cbm = args.perfect_calc_results.cbm
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        try:
            defect_entry = loadfn(d / "defect_entry.json")
            title = defect_entry.name
        except FileNotFoundError:
            title = "No name"
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

        BandEdgeStates(edge_states).to_yaml(d / "band_edge_states.yaml")


def make_defect_formation_energy(args):
    formula = args.perfect_calc_results.structure.composition.reduced_formula
    chem_pot_diag = ChemPotDiag.from_yaml(args.cpd_yaml)
    pcr = args.perfect_calc_results

    defects, defect_entries, corrections, edge_states = [], [], [], []
    for d in args.dirs:
        if args.skip_shallow:
            edge_states = BandEdgeStates.from_yaml(d / "band_edge_states.yaml")
            if edge_states.is_shallow:
                continue
        defects.append(loadfn(d / "calc_results.json"))
        defect_entries.append(loadfn(d / "defect_entry.json"))
        corrections.append(loadfn(d / "correction.json"))

    if args.web_gui:
        from crystal_toolkit.settings import SETTINGS
        import dash_html_components as html
        from crystal_toolkit.helpers.layouts import Column
        import crystal_toolkit.components as ctc
        import dash

        edge_states = []
        for d in args.dirs:
            edge_states.append(BandEdgeStates.from_yaml(d / "band_edge_states.yaml"))

        app = dash.Dash(__name__,
                        suppress_callback_exceptions=True,
                        assets_folder=SETTINGS.ASSETS_PATH,
                        external_stylesheets=['https://codepen.io/chriddyp/pen/bWLwgP.css'])

        cpd_plot_info = CpdPlotInfo(chem_pot_diag)
        cpd_e_component = CpdEnergyComponent(cpd_plot_info,
                                             pcr,
                                             defects,
                                             defect_entries,
                                             corrections,
                                             args.unitcell.vbm,
                                             args.unitcell.cbm,
                                             edge_states)
        my_layout = html.Div([Column(cpd_e_component.layout)])
        ctc.register_crystal_toolkit(app=app, layout=my_layout, cache=None)
        app.run_server(port=args.port)
        return

    abs_chem_pot = chem_pot_diag.abs_chem_pot_dict(args.label)
    title = " ".join([latexify(formula), "point", args.label])
    defect_energies = make_energies(pcr, defects, defect_entries, corrections, abs_chem_pot)

    if args.print:
        defect_energies = slide_energy(defect_energies, args.unitcell.vbm)
        print("         charge          E_f   correction    ")
        for e in defect_energies:
            print(e)
            print("")

        print("-- cross points -- ")
        for e in defect_energies:
            print(e.name)
            print(e.cross_points(ef_min=args.unitcell.vbm,
                                 ef_max=args.unitcell.cbm,
                                 base_ef=args.unitcell.vbm))
            print("")
        return

    plotter = DefectEnergyMplPlotter(title=title,
                                     defect_energies=defect_energies,
                                     vbm=args.unitcell.vbm,
                                     cbm=args.unitcell.cbm,
                                     supercell_vbm=pcr.vbm,
                                     supercell_cbm=pcr.cbm,
                                     y_range=args.y_range,
                                     supercell_edge=args.supercell_edge,
                                     label_line=args.label_line,
                                     add_charges=args.add_charges)

    plotter.construct_plot()
    plotter.plt.savefig(f"energy_{args.label}.pdf")
