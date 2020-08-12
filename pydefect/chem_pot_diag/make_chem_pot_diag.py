# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys
from pathlib import Path
from typing import List, Optional

from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CompositionEnergy, \
    CpdPlotInfo
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiagPlotly3DMplPlotter, \
    ChemPotDiagPlotly2DMplPlotter
from pydefect.util.mp_tools import MpQuery
from pymatgen import Composition, Element
from vise.util.enum import ExtendedEnum

parent = Path(__file__).parent

mp_energies = loadfn(parent / "datasets/mp_atom_energy.yaml")


class AtomEnergyType(ExtendedEnum):
    pbesol = "pbesol"
    hse = "hse"

    @property
    def energies(self):
        if self.pbesol:
            return loadfn(parent / "datasets/vise_pbesol_atom_energy.yaml")


def make_chem_pot_diag_from_mp(elements: List[str],
                               target: Composition,
                               vertex_elements: List[str] = None,
                               atom_energy_yaml: Optional[str] = None):
    """Obtain the energies from Materials Project."""
    properties = ["task_id", "full_formula", "final_energy"]
    query = MpQuery(elements, properties=properties)
    comp_es = set()
    for m in query.materials:
        energy = m["final_energy"]
        if atom_energy_yaml:
            if ".yaml" in atom_energy_yaml:
                energies = loadfn(atom_energy_yaml)
            else:
                energies = AtomEnergyType.from_string(atom_energy_yaml).energies
            diff = {e: energies[e] - mp_energies[e] for e in elements}
            for k, v in Composition(m["full_formula"]).as_dict().items():
                energy += diff[k] * v
        comp_es.add(CompositionEnergy(
            Composition(m["full_formula"]), energy, m["task_id"]))

    return ChemPotDiag(comp_es, target, vertex_elements)


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-e", "--elements", type=str, nargs="+")
    parser.add_argument("-he", "--host_elements", type=Element, nargs="+")
    parser.add_argument("-t", "--target", type=Composition)
    parser.add_argument("-f", "--functional", type=str)
    return parser.parse_args(args)


if __name__ == "__main__":
    x = parse_args(sys.argv[1:])
    cpd = make_chem_pot_diag_from_mp(x.elements, x.target,
                                     x.host_elements, x.functional)
    cpd.to_yaml()
    if cpd.dim == 2:
        plotter = ChemPotDiagPlotly2DMplPlotter(CpdPlotInfo(cpd))
    elif cpd.dim == 3:
        plotter = ChemPotDiagPlotly3DMplPlotter(CpdPlotInfo(cpd))
    plotter.figure.show()
