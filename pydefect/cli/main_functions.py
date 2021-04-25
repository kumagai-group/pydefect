# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from monty.serialization import loadfn
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.defect_energy_plotter import DefectEnergyMplPlotter
from pydefect.analyzer.defect_structure_info import make_defect_structure_info
from pydefect.analyzer.make_band_edge_states import make_band_edge_states
from pydefect.analyzer.make_defect_energy_info import make_defect_energy_info
from pydefect.analyzer.make_defect_energy_summary import \
    make_defect_energy_summary
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    RelativeEnergies, ChemPotDiagMaker, TargetVertices
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiag2DMplPlotter, \
    ChemPotDiag3DMplPlotter
from pydefect.cli.main_tools import sanitize_matrix
from pydefect.cli.vasp.make_efnv_correction import make_efnv_correction
from pydefect.corrections.site_potential_plotter import SitePotentialMplPlotter
from pydefect.input_maker.add_interstitial import append_interstitial
from pydefect.input_maker.defect_set_maker import DefectSetMaker
from pydefect.input_maker.supercell_maker import SupercellMaker
from pymatgen import Composition
from vise.util.logger import get_logger


logger = get_logger(__name__)


def make_standard_and_relative_energies(args):
    comp_energies = CompositionEnergies.from_yaml(
        args.composition_energies_yaml)
    std_energies, rel_energies = comp_energies.std_rel_energies
    std_energies.to_yaml_file()
    rel_energies.to_yaml_file()


def make_cpd_and_vertices(args):
    rel_energies = RelativeEnergies.from_yaml(args.rel_energy_yaml)
    if args.target:
        elements = args.elements or \
                   Composition(args.target).chemical_system.split("-")
    else:
        elements = args.elements or list(rel_energies.all_element_set)

    cpd_maker = ChemPotDiagMaker(rel_energies, elements, args.target)
    cpd = cpd_maker.chem_pot_diag
    cpd.to_json_file()
    cpd.to_target_vertices.to_yaml_file()


def plot_chem_pot_diag(args):
    cpd = args.chem_pot_diag
    if cpd.dim == 2:
        plotter = ChemPotDiag2DMplPlotter(cpd)
    elif cpd.dim == 3:
        plotter = ChemPotDiag3DMplPlotter(cpd)
    else:
        raise ValueError(f"Only 2 or 3 dimensions are supported. "
                         f"Now {cpd.dim} dimensions.")

    plt = plotter.draw_diagram()
    plt.savefig(fname="cpd.pdf")
    plt.show()


def make_supercell(args):
    if args.matrix:
        matrix = sanitize_matrix(args.matrix)
        maker = SupercellMaker(args.unitcell, matrix)
    else:
        kwargs = {}
        if args.min_num_atoms:
            kwargs["min_num_atoms"] = args.min_num_atoms
        if args.max_num_atoms:
            kwargs["max_num_atoms"] = args.max_num_atoms
        maker = SupercellMaker(args.unitcell, **kwargs)

    maker.supercell.structure.to(filename="SPOSCAR")
    maker.supercell_info.to_json_file()


def append_interstitial_to_supercell_info(args):
    supercell_info = append_interstitial(args.supercell_info,
                                         args.base_structure,
                                         args.frac_coords)
    supercell_info.to_json_file()


def pop_interstitial_from_supercell_info(args):
    supercell_info = args.supercell_info
    assert args.index > 0
    supercell_info.interstitials.pop(args.index - 1)
    supercell_info.to_json_file()


def make_defect_set(args):
    supercell_info = loadfn("supercell_info.json")
    if args.oxi_states:
        oxi_states = dict(zip(args.oxi_states[::2], args.oxi_states[1::2]))
    else:
        oxi_states = None
    maker = DefectSetMaker(supercell_info,
                           oxi_states,
                           args.dopants,
                           keywords=args.kwargs)
    maker.defect_set.to_yaml()


def calc_defect_structure_info(args):
    supercell_info = args.supercell_info
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        calc_results = loadfn(d / "calc_results.json")
        defect_entry = loadfn(d / "defect_entry.json")
        defect_str_info = make_defect_structure_info(
            supercell_info.structure,
            defect_entry.structure,
            calc_results.structure,
            dist_tol=args.dist_tolerance,
            symprec=args.symprec,
            init_site_sym=defect_entry.site_symmetry)
        defect_str_info.to_json_file(str(d / "defect_structure_info.json"))


def make_efnv_correction_main_func(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        defect_entry = loadfn(d / "defect_entry.json")
        calc_results = loadfn(d / "calc_results.json")
        efnv = make_efnv_correction(defect_entry.charge,
                                    calc_results,
                                    args.perfect_calc_results,
                                    args.unitcell.dielectric_constant)
        efnv.to_json_file(d / "correction.json")

        title = defect_entry.full_name
        plotter = SitePotentialMplPlotter.from_efnv_corr(
            title=title, efnv_correction=efnv)
        plotter.construct_plot()
        plotter.plt.savefig(fname=d / "correction.pdf")
        plotter.plt.clf()


def make_band_edge_states_main_func(args):
    for d in args.dirs:
        try:
            logger.info(f"Parsing data in {d} ...")
            orb_infos = loadfn(d / "band_edge_orbital_infos.json")
            band_edge_states = make_band_edge_states(orb_infos, args.p_state)
            band_edge_states.to_json_file(d / "band_edge_states.json")
        except ValueError as e:
            print(e)


def make_defect_energy_infos_main_func(args):
    for d in args.dirs:
        logger.info(f"Parsing data in {d} ...")
        defect_entry = loadfn(d / "defect_entry.json")
        calc_results = loadfn(d / "calc_results.json")
        correction = loadfn(d / "correction.json")
        try:
            band_edge_states = loadfn(d / "band_edge_states.json")
        except FileNotFoundError:
            band_edge_states = None
        defect_energy_info = make_defect_energy_info(
            defect_entry=defect_entry,
            calc_results=calc_results,
            correction=correction,
            perfect_calc_results=args.perfect_calc_results,
            standard_energies=args.std_energies,
            band_edge_states=band_edge_states)
        defect_energy_info.to_yaml_file(d / "defect_energy_info.yaml")


def make_defect_energy_summary_main_func(args):
    energy_infos = []
    for d in args.dirs:
        energy_infos.append(DefectEnergyInfo.from_yaml(d / "defect_energy_info.yaml"))
    target_vertices = TargetVertices.from_yaml(args.target_vertices_yaml)
    defect_energy_summary = make_defect_energy_summary(
        energy_infos, target_vertices, args.unitcell, args.p_state)
    defect_energy_summary.to_json_file()


def plot_defect_energy(args):
    plotter = DefectEnergyMplPlotter(
        defect_energy_summary=args.defect_energy_summary,
        chem_pot_label=args.label,
        y_range=args.y_range,
        allow_shallow=args.allow_shallow,
        with_corrections=args.with_corrections,
        label_line=args.label_line,
        add_charges=args.add_charges)
    plotter.construct_plot()
    plotter.plt.savefig(f"energy_{args.label}.pdf")
