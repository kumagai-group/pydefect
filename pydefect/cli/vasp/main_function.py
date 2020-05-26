# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.io.vasp import Vasprun, Outcar
from vise.util.logger import get_logger

from pydefect.cli.main_tools import sanitize_matrix
from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp
from pydefect.cli.vasp.make_unitcell import make_unitcell_from_vasp
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entries_maker import DefectEntriesMaker
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_set_maker import DefectSetMaker
from pydefect.input_maker.supercell_maker import SupercellMaker

logger = get_logger(__name__)


def make_unitcell(args):
    unitcell = make_unitcell_from_vasp(vasprun_band=args.vasprun_band,
                                       outcar_band=args.outcar_band,
                                       outcar_dielectric=args.outcar_dielectric)
    unitcell.to_json_file()


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
    supercell_info = loadfn("supercell_info.json")
    defect_set = DefectSet.from_yaml()
    maker = DefectEntriesMaker(supercell_info, defect_set)
    for defect_entry in maker.defect_entries:
        dir_path = Path(defect_entry.full_name)
        dir_path.mkdir()
        defect_entry.perturbed_structure.to(filename=dir_path / "POSCAR")
        defect_entry.to_json_file(filename=dir_path / "defect_entry.json")


def make_calc_results(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        try:
            dft_results = make_calc_results_from_vasp(
                vasprun=Vasprun(d / defaults.vasprun),
                outcar=Outcar(d / defaults.outcar))
        except IOError:
            logger.warning(f"Parsing data in {d} failed.")
            continue
        dft_results.to_json_file(filename=Path(d) / "calc_results.json")

