# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pydefect import __version__
from pydefect.cli.vasp.main_vasp_util_functions import \
    make_gkfo_correction_from_vasp, \
    calc_defect_charge_info
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args_main_vasp_util(args):
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.
    This command provide some utilities related to the VASP calculations""",
        epilog=f"Author: Yu Kumagai Version: {__version__}")

    subparsers = parser.add_subparsers()

    # ++ parent parser: dirs
    dirs_parser = argparse.ArgumentParser(
        description="", add_help=False)
    dirs_parser.add_argument(
        "-d", "--dirs", nargs="+", type=Path,
        help="Directory paths to be parsed.")

    # ++ parent parser: supercell_info
    si_parser = argparse.ArgumentParser(
        description="", add_help=False)
    si_parser.add_argument(
        "-s", "--supercell_info", required=True, type=loadfn,
        default="supercell_info.json",
        help="Path to the supercell_info.json file.")

    # ++ parent parser: pcr
    pcr_parser = argparse.ArgumentParser(
        description="", add_help=False)
    pcr_parser.add_argument(
        "-pcr", "--perfect_calc_results", required=True, type=loadfn,
        help="Path to the calc_results.json for the perfect supercell.")

    # ++ parent parser: unitcell
    unitcell_parser = argparse.ArgumentParser(
        description="", add_help=False)
    unitcell_parser.add_argument(
        "-u", "--unitcell", required=True, type=loadfn,
        help="Path to the unitcell.json file.")

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

    # -- gkfo correction ------------------------------------------------
    parser_gkfo = subparsers.add_parser(
        name="gkfo",
        description="Generate GKFO correction files.",
        parents=[unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['gkfo'])

    parser_gkfo.add_argument(
        "-iefnv", "--initial_efnv_correction", required=True, type=loadfn,
        help="Path to the initial efnv correction.json file.")
    parser_gkfo.add_argument(
        "-icr", "--initial_calc_results", required=True, type=loadfn,
        help="Path to the initial calc_results.json file.")
    parser_gkfo.add_argument(
        "-fcr", "--final_calc_results", required=True, type=loadfn,
        help="Path to the final calc_results.json file.")
    parser_gkfo.add_argument(
        "-cd", "--charge_diff", required=True, type=int,
        help="Charge difference of final state from initial state.")

    parser_gkfo.set_defaults(func=make_gkfo_correction_from_vasp)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main_vasp_util(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


