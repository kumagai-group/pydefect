# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pydefect import __version__
from pydefect.cli.vasp.main_vasp_functions import make_defect_entries, \
    make_unitcell, make_competing_phase_dirs, \
    make_calc_results, \
    make_band_edge_orb_infos_and_eigval_plot, make_perfect_band_edge_state
from pydefect.defaults import defaults
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args_main_vasp(args):

    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.
    The command is used for creating and analyzing VASP IO files.""",
        epilog=f"Author: Yu Kumagai Version: {__version__}")

    subparsers = parser.add_subparsers()

    # ++ parent parser: dirs
    dirs_parser = argparse.ArgumentParser(
        description="", add_help=False)
    dirs_parser.add_argument(
        "-d", "--dirs", nargs="+", type=Path,
        help="Directory paths to be parsed.")

    # ++ parent parser: pcr
    pcr_parser = argparse.ArgumentParser(
        description="", add_help=False)
    pcr_parser.add_argument(
        "-pcr", "--perfect_calc_results", required=True, type=loadfn,
        help="Path to the calc_results.json for the perfect supercell.")

    # -- unitcell ------------------------------------------------
    parser_unitcell = subparsers.add_parser(
        name="unitcell",
        description="Parse unitcell info and create unitcell.yaml",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['u'])

    parser_unitcell.add_argument(
        "-vb", "--vasprun_band", required=True, type=Vasprun,
        help="vasprun.xml file of band structure calculation.")
    parser_unitcell.add_argument(
        "-ob", "--outcar_band", required=True, type=Outcar,
        help="OUTCAR file of band structure calculation.")
    parser_unitcell.add_argument(
        "-odc", "--outcar_dielectric_clamped", required=True, type=Outcar,
        help="OUTCAR file for ion-clamped dielectric constant.")
    parser_unitcell.add_argument(
        "-odi", "--outcar_dielectric_ionic", required=True, type=Outcar,
        help="OUTCAR file for dielectric constant of ionic contribution.")
    parser_unitcell.add_argument(
        "-n", "--name", type=str, help="System name.")
    parser_unitcell.set_defaults(func=make_unitcell)

    # -- make_poscars ------------------------------------------------
    parser_make_poscars = subparsers.add_parser(
        name="make_poscars",
        description="Make poscar files in directories retrieved from "
                    "the Materials Project.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mp'])

    parser_make_poscars.add_argument(
        "-e", "--elements", required=True, nargs="+", type=str,
        help="Element names considered in chemical potential diagram. ")
    parser_make_poscars.add_argument(
        "--e_above_hull", default=defaults.e_above_hull, type=float,
        help="Energy above hull in eV/atom.")

    parser_make_poscars.set_defaults(func=make_competing_phase_dirs)

    # -- defect_entries ------------------------------------------------
    parser_defect_entries = subparsers.add_parser(
        name="defect_entries",
        description="Make POSCAR files and defect entries in defect "
                    "directories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    parser_defect_entries.set_defaults(func=make_defect_entries)

    # -- calc_results ------------------------------------------------
    parser_calc_results = subparsers.add_parser(
        name="calc_results",
        description="Generate calc_results.json files. "
                    "Run after VASP calculations.",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cr'])
    parser_calc_results.set_defaults(func=make_calc_results)

    # -- band edge orbital infos  ----------------------------------------------
    parser_band_edge_orb_infos = subparsers.add_parser(
        name="band_edge_orbital_infos",
        description="Generate band_edge_orbital_infos.json and "
                    "eigenvalues.pdf.",
        parents=[dirs_parser, pcr_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['beoi'])

    parser_band_edge_orb_infos.set_defaults(
        func=make_band_edge_orb_infos_and_eigval_plot)

    # -- perfect band edge state  ----------------------------------------------
    parser_perf_band_edge_state = subparsers.add_parser(
        name="perfect_band_edge_state",
        description="Generate band_edge_orbital_infos.json and "
                    "eigenvalues.pdf.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pbes'])
    parser_perf_band_edge_state.add_argument(
        "-d", "--dir", type=Path,
        help="Directory path to the perfect supercell calculation.")

    parser_perf_band_edge_state.set_defaults(func=make_perfect_band_edge_state)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main_vasp(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


