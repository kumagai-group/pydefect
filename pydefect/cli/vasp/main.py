# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys

from pymatgen import IStructure

from pydefect.cli.vasp.main_function import make_supercell
from pydefect.version import __version__


def parse_args(args):

    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.""",
        epilog=f"""                                 
    Author: Yu Kumagai
    Version: {__version__}""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers()

    # -- make_supercell ------------------------------------------------
    parser_make_supercell = subparsers.add_parser(
        name="make_supercell",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ms'])

    parser_make_supercell.add_argument(
        "-p", "--unitcell",
        default="POSCAR",
        type=IStructure.from_file)
    parser_make_supercell.add_argument(
        "--matrix",
        default=None,
        nargs="+",
        type=int)
    parser_make_supercell.add_argument(
        "--min_num_atoms",
        default=50,
        type=int)
    parser_make_supercell.add_argument(
        "--max_num_atoms",
        default=300,
        type=int)

    parser_make_supercell.set_defaults(func=make_supercell)

    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()

