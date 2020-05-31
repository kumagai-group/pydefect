# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
from pathlib import Path

from monty.serialization import loadfn
from pymatgen import IStructure, Composition
from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.main_function import make_supercell, make_defect_set, \
    make_defect_entries, make_unitcell, make_competing_phase_dirs, \
    make_chem_pot_diag, make_efnv_correction_from_vasp, print_file, \
    make_calc_results, make_defect_formation_energy
from pydefect.defaults import defaults
from pydefect.version import __version__


def parse_args(args):

    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.""",
        epilog=f"Author: Yu Kumagai Version: {__version__}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    # ++ parent parser: dirs
    dirs_parser = argparse.ArgumentParser(
        description="", add_help=False)
    dirs_parser.add_argument(
        "-d", "--dirs",
        required=True,
        nargs="+",
        type=Path)

    # ++ parent parser: refs
    refs_parser = argparse.ArgumentParser(
        description="", add_help=False)
    refs_parser.add_argument(
        "-pcr", "--perfect_calc_results",
        required=True,
        type=loadfn)
    refs_parser.add_argument(
        "-u", "--unitcell",
        required=True,
        type=loadfn)

    # -- print ------------------------------------------------
    parser_print = subparsers.add_parser(
        name="print",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['p'])

    parser_print.add_argument(
        "-f", "--filename",
        dest="obj",
        required=True,
        type=loadfn)

    parser_print.set_defaults(func=print_file)

    # -- unitcell ------------------------------------------------
    parser_unitcell = subparsers.add_parser(
        name="unitcell",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['u'])

    parser_unitcell.add_argument(
        "-vb", "--vasprun_band",
        required=True,
        type=Vasprun)
    parser_unitcell.add_argument(
        "-ob", "--outcar_band",
        required=True,
        type=Outcar)
    parser_unitcell.add_argument(
        "-od", "--outcar_dielectric",
        required=True,
        type=Outcar)

    parser_unitcell.set_defaults(func=make_unitcell)

    # -- make_poscars ------------------------------------------------
    parser_make_poscars = subparsers.add_parser(
        name="make_poscars",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mp'])

    parser_make_poscars.add_argument(
        "-e", "--elements",
        required=True,
        nargs="+",
        type=str)
    parser_make_poscars.add_argument(
        "--e_above_hull",
        default=defaults.e_above_hull,
        type=float)

    parser_make_poscars.set_defaults(func=make_competing_phase_dirs)

    # -- cpd ------------------------------------------------
    parser_cpd = subparsers.add_parser(
        name="cpd",
        description="",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cpd'])

    parser_cpd.add_argument(
        "-t", "--target",
        type=Composition)

    parser_cpd.set_defaults(func=make_chem_pot_diag)

    # -- supercell ------------------------------------------------
    parser_supercell = subparsers.add_parser(
        name="supercell",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['s'])

    parser_supercell.add_argument(
        "-p", "--unitcell",
        default="POSCAR",
        type=IStructure.from_file)
    parser_supercell.add_argument(
        "--matrix",
        default=None,
        nargs="+",
        type=int)
    parser_supercell.add_argument(
        "--min_num_atoms",
        default=50,
        type=int)
    parser_supercell.add_argument(
        "--max_num_atoms",
        default=300,
        type=int)

    parser_supercell.set_defaults(func=make_supercell)

    # -- defect_set ------------------------------------------------
    parser_defect_set = subparsers.add_parser(
        name="defect_set",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ds'])

    parser_defect_set.add_argument(
        "-o", "--oxi_states",
        nargs="+",
        type=str)
    parser_defect_set.add_argument(
        "-d", "--dopants",
        nargs="+",
        type=str)
    parser_defect_set.add_argument(
        "-k", "--kwargs",
        nargs="+",
        type=str)

    parser_defect_set.set_defaults(func=make_defect_set)

    # -- defect_entries ------------------------------------------------
    parser_defect_entries = subparsers.add_parser(
        name="defect_entries",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    parser_defect_entries.set_defaults(func=make_defect_entries)

    # -- calc_results ------------------------------------------------
    parser_calc_results = subparsers.add_parser(
        name="calc_results",
        description="",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cr'])
    parser_calc_results.set_defaults(func=make_calc_results)

    # -- efnv correction ------------------------------------------------
    parser_efnv = subparsers.add_parser(
        name="efnv",
        description="",
        parents=[dirs_parser, refs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['efnv'])

    parser_efnv.set_defaults(func=make_efnv_correction_from_vasp)

    # -- defect formation energy ----------------------------------------------
    parser_efnv = subparsers.add_parser(
        name="defect_formation_energy",
        description="",
        parents=[dirs_parser, refs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['e'])

    parser_efnv.add_argument(
        "-c", "--chem_pot_diag",
        required=True,
        type=loadfn)
    parser_efnv.add_argument(
        "-l", "--label",
        required=True,
        type=str)
    parser_efnv.add_argument(
        "-y", "--y_range",
        nargs=2,
        type=float)

    parser_efnv.set_defaults(func=make_defect_formation_energy)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


