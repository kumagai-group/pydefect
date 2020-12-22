# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.


import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pydefect import __version__
from pydefect.cli.main_tools import str_int_to_int
from pydefect.cli.vasp.main_function import make_supercell, make_defect_set, \
    make_defect_entries, make_unitcell, make_competing_phase_dirs, \
    make_chem_pot_diag, make_efnv_correction_from_vasp, print_file, \
    make_calc_results, make_defect_formation_energy, make_defect_eigenvalues, \
    make_edge_characters, make_edge_states, \
    append_interstitial_to_supercell_info, pop_interstitial_from_supercell_info, \
    plot_chem_pot_diag, make_gkfo_correction_from_vasp
from pydefect.defaults import defaults
from pymatgen import IStructure, Composition, Structure
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)


def parse_args(args):

    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package that helps researchers to do first-principles point 
    defect calculations with the VASP code.""",
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

    # ++ parent parser: unitcell
    unitcell_parser = argparse.ArgumentParser(
        description="", add_help=False)
    unitcell_parser.add_argument(
        "-u", "--unitcell", required=True, type=loadfn,
        help="Path to the unitcell.json file.")

    # -- print ------------------------------------------------
    parser_print = subparsers.add_parser(
        name="print",
        description="Print the json file info.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['p'])

    parser_print.add_argument(
        "-f", "--filename", dest="obj", required=True, type=loadfn,
        help="Path to the json file to be printed.")

    parser_print.set_defaults(func=print_file)

    # -- unitcell ------------------------------------------------
    parser_unitcell = subparsers.add_parser(
        name="unitcell",
        description="Parse unitcell info and create unitcell.json",
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

    # -- make_cpd ------------------------------------------------
    parser_mcpd = subparsers.add_parser(
        name="make_cpd",
        parents=[dirs_parser],
        description="Make chemical potential diagram.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcpd'])

    parser_mcpd.add_argument(
        "-y", "--yaml", type=str, default="cpd.yaml",
        help="Saved yaml file name.")
    parser_mcpd.add_argument(
        "-t", "--target", type=Composition, required=True,
        help="Target composition, e.g., MgO.")
    parser_mcpd.add_argument(
        "-e", "--elements", type=str, nargs="+",
        help="Element names considered in chemical potential diagram. Used for"
             "creating the diagram from mp data with atom energy alignment.")
    parser_mcpd.add_argument(
        "-a", "--atom_energy_yaml", type=str,
        help="Yaml file storing atom energies for energy alignment.")
    parser_mcpd.add_argument(
        "-u", "--update", action="store_true",
        help="Used when one wants to update the cpd.yaml with calculation "
             "results.")

    parser_mcpd.set_defaults(func=make_chem_pot_diag)

    # -- plot_cpd ------------------------------------------------
    parser_pcpd = subparsers.add_parser(
        name="plot_cpd",
        description="Plot chemical potential diagram and make pdf file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pcpd'])

    parser_pcpd.add_argument(
        "-y", "--yaml", default="cpd.yaml", type=str,
        help="cpd.yaml file name")

    parser_pcpd.set_defaults(func=plot_chem_pot_diag)

    # -- supercell ------------------------------------------------
    parser_supercell = subparsers.add_parser(
        name="supercell",
        description="Make SPOSCAR file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['s'])

    parser_supercell.add_argument(
        "-p", "--unitcell", default="POSCAR", type=IStructure.from_file,
        help="Base POSCAR, which must be the standardized primitive cell.")
    parser_supercell.add_argument(
        "--matrix", nargs="+", type=int,
        help="Supercell matrix applied to the conventional cell. "
             "1, 3 or 9 components are accepted.")
    parser_supercell.add_argument(
        "--min_atoms", dest="min_num_atoms", default=50, type=int,
        help="Minimum number of atoms")
    parser_supercell.add_argument(
        "--max_atoms", dest="max_num_atoms", default=300, type=int,
        help="Maximum number of atoms")

    parser_supercell.set_defaults(func=make_supercell)

    # -- defect_set ------------------------------------------------
    parser_defect_set = subparsers.add_parser(
        name="defect_set",
        description="Make defect_in.yaml file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ds'])

    parser_defect_set.add_argument(
        "-o", "--oxi_states", nargs="+", type=str_int_to_int,
        help="Oxidation states in integers, e.g., Mg 2 O -2.")
    parser_defect_set.add_argument(
        "-d", "--dopants", nargs="+", type=str,
        help="Names of dopant element, e.g., Al Ga.")
    parser_defect_set.add_argument(
        "-k", "--kwargs", nargs="+", type=str,
        help="Keywords used to screen the target defects. "
             "Regular expression can be used.")

    parser_defect_set.set_defaults(func=make_defect_set)

    # -- defect_entries ------------------------------------------------
    parser_defect_entries = subparsers.add_parser(
        name="defect_entries",
        description="Make POSCAR files and defect entries in defect "
                    "directories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['de'])

    parser_defect_entries.set_defaults(func=make_defect_entries)

    # -- append_interstitial ------------------------------------------------
    parser_append_interstitial = subparsers.add_parser(
        name="append_interstitial",
        description="Append interstitial information to supercell_info.yaml",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ai'])

    parser_append_interstitial.add_argument(
        "-s", "--supercell_info", required=True, type=loadfn,
        help="Path to the supercell_info.json file.")
    parser_append_interstitial.add_argument(
        "-p", "--base_structure", required=True, type=Structure.from_file,
        help="POSCAR file used for generating AECCAR files, which must be the "
             "standardized primitive cell.")
    parser_append_interstitial.add_argument(
        "-c", "--frac_coords", required=True, nargs=3, type=float,
        help="Fractional coordinates in the standardized primitive cell added "
             "as an interstitial ")

    parser_append_interstitial.set_defaults(
        func=append_interstitial_to_supercell_info)

    # -- pop_interstitial ------------------------------------------------
    parser_pop_interstitial = subparsers.add_parser(
        name="pop",
        description="Pop interstitial site from supercell_info.json. ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pi'])

    parser_pop_interstitial.add_argument(
        "-s", "--supercell_info", required=True, type=loadfn,
        help="Path to the supercell_info.json file.")
    parser_pop_interstitial.add_argument(
        "-i", "--index", required=True, type=int,
        help="Popped interstitial index. Note that the indices begin from 1.")

    parser_pop_interstitial.set_defaults(
        func=pop_interstitial_from_supercell_info)

    # -- calc_results ------------------------------------------------
    parser_calc_results = subparsers.add_parser(
        name="calc_results",
        description="Generate calc_results.json files. "
                    "Run after VASP calculations.",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cr'])
    parser_calc_results.set_defaults(func=make_calc_results)

    # -- efnv correction ------------------------------------------------
    parser_efnv = subparsers.add_parser(
        name="efnv",
        description="Generate extended FNV-correction related efnv_correction "
                    "files. ",
        parents=[dirs_parser, pcr_parser, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['efnv'])

    parser_efnv.set_defaults(func=make_efnv_correction_from_vasp)

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

    # -- defect eigenvalues ------------------------------------------------
    parser_eig = subparsers.add_parser(
        name="eig",
        description="Generate plot of eigenvalues and occupations.",
        parents=[dirs_parser, pcr_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['eig'])

    parser_eig.set_defaults(func=make_defect_eigenvalues)

    # -- band edge characters ------------------------------------------------
    parser_ec = subparsers.add_parser(
        name="make_edge_characters",
        description="Make edge_characters.json files that are used for "
                    "analyzing the eigenvalues.",
        parents=[dirs_parser, pcr_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mec'])

    parser_ec.set_defaults(func=make_edge_characters)

    # -- band edge states ------------------------------------------------
    parser_es = subparsers.add_parser(
        name="edge_states",
        description="Show edge state for each spin channel.",
        parents=[dirs_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['es'])

    parser_es.add_argument(
        "-p", "--perfect_edge_characters", required=True, type=loadfn,
        help="Path to the edge_characters.json for the perfect supercell.")

    parser_es.set_defaults(func=make_edge_states)

    # -- defect formation energy ----------------------------------------------
    parser_energy = subparsers.add_parser(
        name="defect_formation_energy",
        description="Show and plot defect formation energies.",
        parents=[dirs_parser, pcr_parser, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['e'])

    parser_energy.add_argument(
        "-c", "--cpd_yaml", required=True, type=str,
        help="cpd.yaml file path.")
    parser_energy.add_argument(
        "-l", "--label", type=str, help="Label in cpd.pdf")
    parser_energy.add_argument(
        "-y", "--y_range", nargs=2, type=float,
        help="Energy range in y-axis")
    parser_energy.add_argument(
        "-s", "--skip_shallow", action="store_true",
        help="Set when the energies of shallow defects are omitted.")
    parser_energy.add_argument(
        "-p", "--print", action="store_true",
        help="Print the defect formation energy as a list.")
    parser_energy.add_argument(
        "--supercell_edge", action="store_true",
        help="Whether set the band edges to supercell ones instead of "
             "unitcell ones.")
    parser_energy.add_argument(
        "--no_label_line", dest="label_line", action="store_false",
        help="Set the labels not to locate on the lines.")
    parser_energy.add_argument(
        "--no_add_charges", dest="add_charges", action="store_false",
        help="Set the charges not to appear.")
    parser_energy.add_argument(
        "-wg", "--web_gui", action="store_true",
        help="Show web GUI.")
    parser_energy.add_argument(
        "--port", type=int, default=8050,
        help="Port for web GUI.")

    parser_energy.set_defaults(func=make_defect_formation_energy)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()


