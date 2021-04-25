# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
import warnings
from pathlib import Path

from pydefect.cli.main import epilog, description
from pydefect.cli.vasp.main_vasp_util_functions import \
    calc_defect_charge_info
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args_main_vasp_util(args):
    parser = argparse.ArgumentParser(epilog=epilog,
                                     description=description + """      
    This command provides some utilities related to the VASP calculations""")

    subparsers = parser.add_subparsers()

    # -- calc defect charge info -----------------------------------------------
    parser_calc_def_charge_info = subparsers.add_parser(
        name="calc_defect_charge_info",
        description="Calc defect charge info.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cdc'])

    parser_calc_def_charge_info.add_argument(
        "-p", "--parchgs", type=str, nargs="+", required=True,
        help="PARCHG files.")
    parser_calc_def_charge_info.add_argument(
        "-v", "--vesta_file", type=Path, help="defect.vesta file")
    parser_calc_def_charge_info.add_argument(
        "-b", "--bin_interval", type=float, default=0.2)
    parser_calc_def_charge_info.add_argument(
        "-g", "--grids_dirname", type=Path)

    parser_calc_def_charge_info.set_defaults(func=calc_defect_charge_info)

    return parser.parse_args(args)


def main():
    args = parse_args_main_vasp_util(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


