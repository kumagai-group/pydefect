# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import argparse
import sys
import warnings

from monty.serialization import loadfn
from pydefect.cli.main import description, epilog, add_sub_parser
from pydefect.cli.main_util_functions import make_gkfo_correction_from_vasp, \
    composition_energies_from_mp, add_interstitials_from_local_extrema, \
    make_defect_vesta_file, show_u_values
from pydefect.defaults import defaults
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args_main_util(args):
    parser = argparse.ArgumentParser(epilog=epilog,
                                     description=description + """                            
    This command provide some utilities related to the VASP calculations""")

    subparsers = parser.add_subparsers()
    dirs_parser = add_sub_parser(argparse, name="dirs")
    unitcell_parser = add_sub_parser(argparse, name="unitcell")
    si_parser = add_sub_parser(argparse, name="supercell_info")

    # -- composition energies from mp ------------------------------------------
    parser_comp_es_from_mp = subparsers.add_parser(
        name="composition_energies_from_mp",
        description="Retrieve composition energies from Materials Project.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cefm'])

    parser_comp_es_from_mp.add_argument(
        "-e", "--elements", type=str, nargs="+",
        help="Element names considered in chemical potential diagram. Used for "
             "creating the diagram from mp data with atom energy alignment.")
    parser_comp_es_from_mp.add_argument(
        "-a", "--atom_energy_yaml", type=str,
        help="Yaml file storing atom energies for energy alignment.")

    parser_comp_es_from_mp.set_defaults(func=composition_energies_from_mp)

    # -- show u values ------------------------------------------
    parser_show_u_values = subparsers.add_parser(
        name="show_u_values",
        description="Show U values for defects.",
        parents=[add_sub_parser(argparse, "defect_energy_summary")],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['u'])

    parser_show_u_values.set_defaults(func=show_u_values)

    # -- add interstitials from local extrema ----------------------------------
    parser_ai = subparsers.add_parser(
        name="add_interstitials_from_local_extrema",
        description="Add interstitials via volumetric_data_local_extrema.json "
                    "file.",
        parents=[si_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ai'])
    parser_ai.add_argument(
        "--local_extrema", required=True, type=loadfn,
        help="volumetric_data_local_extrema.json file name.")
    parser_ai.add_argument(
        "-i", "--indices", required=True, type=int, nargs="+",
        help="Indices starting from 1 to be added to SupercellInfo.")
    parser_ai.set_defaults(func=add_interstitials_from_local_extrema)

    # -- make defect vesta file ------------------------------------------------
    parser_dvf = subparsers.add_parser(
        name="defect_vesta_file",
        description="Generate VESTA file for visualizing defect structure.",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['dvf'])
    parser_dvf.add_argument(
        "--cutoff", type=float, default=defaults.show_structure_cutoff,
        help="Cutoff radius in which atoms are shown in VESTA files.")
    parser_dvf.add_argument(
        "--min_displace_w_arrows", type=float, default=0.1,
        help="Minimum displacement that shows arrows to visualize "
             "displacements.")
    parser_dvf.add_argument(
        "--arrow_factor", type=float, default=3.0,
        help="Factor to enhance the displacements.")
    parser_dvf.add_argument(
        "--title", type=str, help="Title to be shown in VESTA files.")

    parser_dvf.set_defaults(func=make_defect_vesta_file)

    # -- gkfo correction -------------------------------------------------------
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
    args = parse_args_main_util(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


