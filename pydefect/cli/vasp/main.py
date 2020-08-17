# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
from pathlib import Path

from monty.serialization import loadfn
from pydefect.cli.main_tools import str_int_to_int
from pydefect.cli.vasp.main_function import make_supercell, make_defect_set, \
    make_defect_entries, make_unitcell, make_competing_phase_dirs, \
    make_chem_pot_diag, make_efnv_correction_from_vasp, print_file, \
    make_calc_results, make_defect_formation_energy, make_defect_eigenvalues, \
    make_edge_characters, make_edge_states, \
    append_interstitial_to_supercell_info, pop_interstitial_from_supercell_info, \
    plot_chem_pot_diag
from pydefect.defaults import defaults
from pydefect.version import __version__
from pymatgen import IStructure, Composition, Structure
from pymatgen.io.vasp import Vasprun, Outcar


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
        "-d", "--dirs", required=True, nargs="+", type=Path)

    # ++ parent parser: pcr
    pcr_parser = argparse.ArgumentParser(
        description="", add_help=False)
    pcr_parser.add_argument(
        "-pcr", "--perfect_calc_results", required=True, type=loadfn)

    # ++ parent parser: unitcell
    unitcell_parser = argparse.ArgumentParser(
        description="", add_help=False)
    unitcell_parser.add_argument(
        "-u", "--unitcell", required=True, type=loadfn)

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
        description="Parse unitcell info and create unitcell.json",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['u'])

    parser_unitcell.add_argument(
        "-vb", "--vasprun_band", required=True, type=Vasprun)
    parser_unitcell.add_argument(
        "-ob", "--outcar_band", required=True, type=Outcar)
    parser_unitcell.add_argument(
        "-od", "--outcar_dielectric", required=True, type=Outcar)

    parser_unitcell.set_defaults(func=make_unitcell)

    # -- make_poscars ------------------------------------------------
    parser_make_poscars = subparsers.add_parser(
        name="make_poscars",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mp'])

    parser_make_poscars.add_argument(
        "-e", "--elements", required=True, nargs="+", type=str)
    parser_make_poscars.add_argument(
        "--e_above_hull", default=defaults.e_above_hull, type=float)

    parser_make_poscars.set_defaults(func=make_competing_phase_dirs)

    # -- make_cpd ------------------------------------------------
    parser_mcpd = subparsers.add_parser(
        name="make_cpd",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcpd'])

    parser_mcpd.add_argument("-y", "--yaml", type=str, default="cpd.yaml")
    parser_mcpd.add_argument("-t", "--target", type=Composition, required=True)
    parser_mcpd.add_argument("-e", "--elements", type=str, nargs="+")
    parser_mcpd.add_argument("-f", "--functional", type=str)
    parser_mcpd.add_argument("-d", "--dirs", nargs="+", type=Path)
    parser_mcpd.add_argument("-u", "--update", action="store_true")

    parser_mcpd.set_defaults(func=make_chem_pot_diag)

    # -- plot_cpd ------------------------------------------------
    parser_pcpd = subparsers.add_parser(
        name="plot_cpd",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pcpd'])

    parser_pcpd.add_argument("-y", "--yaml", type=str)

    parser_pcpd.set_defaults(func=plot_chem_pot_diag)

    # -- supercell ------------------------------------------------
    parser_supercell = subparsers.add_parser(
        name="supercell",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['s'])

    parser_supercell.add_argument(
        "-p", "--unitcell", default="POSCAR", type=IStructure.from_file)
    parser_supercell.add_argument("--matrix", default=None, nargs="+", type=int)
    parser_supercell.add_argument("--min_num_atoms", default=50, type=int)
    parser_supercell.add_argument("--max_num_atoms", default=300, type=int)

    parser_supercell.set_defaults(func=make_supercell)

    # -- defect_set ------------------------------------------------
    parser_defect_set = subparsers.add_parser(
        name="defect_set",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ds'])

    parser_defect_set.add_argument(
        "-o", "--oxi_states", nargs="+", type=str_int_to_int)
    parser_defect_set.add_argument("-d", "--dopants", nargs="+", type=str)
    parser_defect_set.add_argument("-k", "--kwargs", nargs="+", type=str)

    parser_defect_set.set_defaults(func=make_defect_set)

    # -- defect_entries ------------------------------------------------
    parser_defect_entries = subparsers.add_parser(
        name="defect_entries",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    parser_defect_entries.set_defaults(func=make_defect_entries)

    # -- append_interstitial ------------------------------------------------
    parser_append_interstitial = subparsers.add_parser(
        name="append_interstitial",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ai'])

    parser_append_interstitial.add_argument(
        "-s", "--supercell_info", required=True, type=loadfn)
    parser_append_interstitial.add_argument(
        "-p", "--base_structure", required=True, type=Structure.from_file)
    parser_append_interstitial.add_argument(
        "-c", "--frac_coords", required=True, nargs=3, type=float)

    parser_append_interstitial.set_defaults(func=append_interstitial_to_supercell_info)

    # -- pop_interstitial ------------------------------------------------
    parser_pop_interstitial = subparsers.add_parser(
        name="pop",
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pi'])

    parser_pop_interstitial.add_argument(
        "-s", "--supercell_info", required=True, type=loadfn)
    parser_pop_interstitial.add_argument(
        "-i", "--index", required=True, type=int)

    parser_pop_interstitial.set_defaults(func=pop_interstitial_from_supercell_info)

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
        parents=[dirs_parser, pcr_parser, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['efnv'])

    parser_efnv.set_defaults(func=make_efnv_correction_from_vasp)

    # -- defect eigenvalues ------------------------------------------------
    parser_eig = subparsers.add_parser(
        name="eig",
        description="",
        parents=[dirs_parser, pcr_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['eig'])

    parser_eig.set_defaults(func=make_defect_eigenvalues)

    # -- band edge characters ------------------------------------------------
    parser_ec = subparsers.add_parser(
        name="edge_characters",
        description="",
        parents=[dirs_parser, pcr_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ec'])

    parser_ec.set_defaults(func=make_edge_characters)

    # -- band edge states ------------------------------------------------
    parser_es = subparsers.add_parser(
        name="edge_states",
        description="",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['es'])

    parser_es.add_argument(
        "-p", "--perfect_edge_characters", required=True, type=loadfn)

    parser_es.set_defaults(func=make_edge_states)
    # -- defect formation energy ----------------------------------------------
    parser_energy = subparsers.add_parser(
        name="defect_formation_energy",
        description="",
        parents=[dirs_parser, pcr_parser, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['e'])

    parser_energy.add_argument(
        "-c", "--chem_pot_diag", required=True, type=str)
    parser_energy.add_argument("-l", "--label", required=True, type=str)
    parser_energy.add_argument("-y", "--y_range", nargs=2, type=float)
    parser_energy.add_argument("-s", "--skip_shallow", action="store_true")
    parser_energy.add_argument("-p", "--print", action="store_true")

    parser_energy.set_defaults(func=make_defect_formation_energy)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


