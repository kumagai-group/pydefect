# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    RelativeEnergies, ChemPotDiagMaker, ChemPotDiag
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiag2DPlotlyPlotter, \
    ChemPotDiag3DPlotlyPlotter, ChemPotDiag2DMplPlotter, ChemPotDiag3DMplPlotter
from pymatgen import Composition
from vise.util.logger import get_logger


logger = get_logger(__name__)


def make_std_and_rel_energies(args):
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
        elements = list(rel_energies.all_element_set)

    cpd_maker = ChemPotDiagMaker(rel_energies, elements, args.target)
    cpd = cpd_maker.chem_pot_diag
    cpd.to_json()
    cpd.to_target_vertices.to_yaml_file()


def plot_chem_pot_diag(args):
    cpd: ChemPotDiag = args.chem_pot_diag
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

