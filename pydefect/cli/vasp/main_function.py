# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.io.vasp import Vasprun, Outcar
from vise.util.logger import get_logger

from pydefect.analyzer.defect_energy import make_defect_energies
from pydefect.analyzer.defect_energy_plotter import DefectEnergyPlotter
from pydefect.analyzer.make_defect_energy import make_single_defect_energy
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiag2DPlotter, \
    ChemPotDiag3DPlotter
from pydefect.cli.main_tools import sanitize_matrix
from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp
from pydefect.cli.vasp.make_poscars_from_query import make_poscars_from_query
from pydefect.cli.vasp.make_unitcell import make_unitcell_from_vasp
from pydefect.corrections.efnv_correction.make_efnv_correction import \
    make_efnv_correction
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entries_maker import DefectEntriesMaker
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_set_maker import DefectSetMaker
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.input_maker.supercell_maker import SupercellMaker
from pydefect.util.error_classes import CpdNotSupportedError
from pydefect.util.mp_tools import MpQuery

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
    energies = {}
    for d in args.dirs:
        vasprun = Vasprun(d / defaults.vasprun)
        composition = vasprun.final_structure.composition
        energy = vasprun.final_energy
        energies[str(composition)] = energy

    cpd = ChemPotDiag(energies, args.target)
    cpd.to_json_file()

    if cpd.dim == 2:
        plotter = ChemPotDiag2DPlotter(CpdPlotInfo(cpd))
    elif cpd.dim == 3:
        plotter = ChemPotDiag3DPlotter(CpdPlotInfo(cpd))
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


def make_defect_set(args):
    supercell_info = loadfn("supercell_info.json")
    maker = DefectSetMaker(supercell_info,
                           args.oxi_states,
                           args.dopants,
                           keywords=args.kwargs)
    maker.defect_set.to_yaml()


def make_defect_entries(args):
    supercell_info: SupercellInfo = loadfn("supercell_info.json")
    perfect = Path("perfect")
    perfect.mkdir()

    supercell_info.structure.to(filename=perfect / "POSCAR")
    defect_set = DefectSet.from_yaml()
    logger.info("Making perfect dir...")
    maker = DefectEntriesMaker(supercell_info, defect_set)

    for defect_entry in maker.defect_entries:
        dir_path = Path(defect_entry.full_name)
        logger.info(f"Making {dir_path} dir...")
        dir_path.mkdir()
        defect_entry.perturbed_structure.to(filename=dir_path / "POSCAR")
        defect_entry.to_json_file(filename=dir_path / "defect_entry.json")
        defect_entry.to_prior_info(filename=dir_path / "prior_info.yaml")


def make_calc_results(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        calc_results = make_calc_results_from_vasp(
                vasprun=Vasprun(d / defaults.vasprun),
                outcar=Outcar(d / defaults.outcar))
        calc_results.to_json_file(filename=Path(d) / "calc_results.json")


def make_efnv_correction_from_vasp(args):
    for d in args.dirs:
        defect_entry = loadfn(d / "defect_entry.json")
        calc_results = loadfn(d / "calc_results.json")
        efnv = make_efnv_correction(defect_entry.charge,
                                    calc_results,
                                    args.perfect_calc_results,
                                    args.unitcell.dielectric_constant)
        efnv.to_json_file(d / "correction.json")


def make_defect_formation_energy(args):
    title = args.perfect_calc_results.structure.composition.reduced_formula
    abs_chem_pot = args.chem_pot_diag.abs_chem_pot_dict(args.label)

    single_energies = []
    for d in args.dirs:
        single_energies.append(
            make_single_defect_energy(args.perfect_calc_results,
                                      loadfn(d / "calc_results.json"),
                                      loadfn(d / "defect_entry.json"),
                                      abs_chem_pot))
#                                      loadfn(d / "correction.json")))
    DefectEnergyPlotter(title=title,
                        defect_energies=make_defect_energies(single_energies),
                        vbm=args.unitcell.vbm,
                        cbm=args.unitcell.cbm,
                        supercell_vbm=args.perfect_calc_results.vbm,
                        supercell_cbm=args.perfect_calc_results.cbm,
                        y_range=args.y_range)


