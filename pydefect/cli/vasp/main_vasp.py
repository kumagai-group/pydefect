# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
import warnings
from pathlib import Path

from pydefect.cli.main import epilog, description, add_sub_parser
from pydefect.cli.vasp.main_vasp_functions import make_defect_entries, \
    make_unitcell, make_competing_phase_dirs, \
    make_calc_results, \
    make_band_edge_orb_infos_and_eigval_plot, make_perfect_band_edge_state, \
    make_local_extrema
from pydefect.defaults import defaults
from pymatgen.io.vasp import Vasprun, Outcar, Chgcar
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args_main_vasp(args):

    parser = argparse.ArgumentParser(epilog=epilog,
                                     description=description + """       
    The command is used for creating and analyzing VASP IO files.""")

    subparsers = parser.add_subparsers()
    dirs_parser = add_sub_parser(argparse, name="dirs")
    pbes_parser = add_sub_parser(argparse, name="perfect_band_edge_state")

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

    # -- make_local_extrema ----------------------------------------------------
    parser_make_local_extrema = subparsers.add_parser(
        name="local_extrema",
        description="Make local volumetric_data_local_extrema.json file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['le'])

    parser_make_local_extrema.add_argument(
        "-v", "--volumetric_data", type=Chgcar.from_file, required=True,
        nargs="+",
        help="File names such as CHGCAR or LOCPOT. When multiple files are "
             "provided, the summed data will be parsed.")
    parser_make_local_extrema.add_argument(
        "--disable_find_min", dest="find_min", action="store_false",
        help="Set when local maxima are searched instead of local minima.")
    parser_make_local_extrema.add_argument(
        "-i", "--info", type=str,
        help="Information related to the parsed volumetric data.")
    parser_make_local_extrema.add_argument(
        "--threshold_frac", default=None, type=float,
        help="""Optional fraction of extrema shown, which returns 
        `threshold_frac * tot_num_extrema` extrema fractional coordinates 
         based on highest/lowest intensity. It takes from 0 to 1.,
         If not set, 1 will be set.""")
    parser_make_local_extrema.add_argument(
        "--threshold_abs", default=None, type=float,
        help="""Optional filter. When searching for local
        minima, intensity <= threshold_abs returns; when searching for
        local maxima, intensity >= threshold_abs returns.

        Note that threshold_abs and threshold_frac should not set in the
        same time.""")
    parser_make_local_extrema.add_argument(
        "--min_dist", type=float, metavar="Angstrom", default=0.5,
        help="""Used to remove the predicted sites that are too close to 
        *existing atoms* in the structure. The minimum distance that a vertex 
        needs to be from existing atoms. Set 0 when switch off this flag.""")
    parser_make_local_extrema.add_argument(
        "--tol", type=float, metavar="Angstrom", default=0.5,
        help="""Group interstitials that are too close together using a tol.
        Set 0 when switch off this flag.""")
    parser_make_local_extrema.add_argument(
        "--radius", type=float, default=0.4,
        help="Radius of sphere around each site to evaluate the average "
             "quantity.")

    parser_make_local_extrema.set_defaults(func=make_local_extrema)
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

    # -- band edge orbital infos  ----------------------------------------------
    parser_band_edge_orb_infos = subparsers.add_parser(
        name="band_edge_orbital_infos",
        description="Generate band_edge_orbital_infos.json and "
                    "eigenvalues.pdf.",
        parents=[dirs_parser, pbes_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['beoi'])

    parser_band_edge_orb_infos.set_defaults(
        func=make_band_edge_orb_infos_and_eigval_plot)

    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main_vasp(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


