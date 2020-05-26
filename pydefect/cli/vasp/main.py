# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys

from pymatgen import IStructure
from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.main_function import make_supercell, make_defect_set, \
    make_defect_entries, make_unitcell
from pydefect.version import __version__


def parse_args(args):

    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.""",
        epilog=f"Author: Yu Kumagai Version: {__version__}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

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

    # -- make_defect_set ------------------------------------------------
    parser_make_defect_set = subparsers.add_parser(
        name="defect_set",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ds'])

    parser_make_defect_set.add_argument(
        "-o", "--oxi_states",
        nargs="+",
        type=str)
    parser_make_defect_set.add_argument(
        "-d", "--dopants",
        nargs="+",
        type=str)
    parser_make_defect_set.add_argument(
        "-k", "--kwargs",
        nargs="+",
        type=str)

    parser_make_defect_set.set_defaults(func=make_defect_set)

    # -- make_defect_entries ------------------------------------------------
    parser_make_defect_entries = subparsers.add_parser(
        name="defect_entries",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    parser_make_defect_entries.set_defaults(func=make_defect_entries)
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

