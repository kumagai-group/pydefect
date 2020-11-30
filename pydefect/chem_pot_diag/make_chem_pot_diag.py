# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys
from itertools import groupby
from pathlib import Path
from typing import List, Optional, Union

from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CompositionEnergy, \
    CpdPlotInfo
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiagPlotly3DMplPlotter, \
    ChemPotDiagPlotly2DMplPlotter
from pydefect.util.mp_tools import MpQuery
from pymatgen import Composition, Element
from vise.atom_energies.atom_energy import mp_energies
from vise.util.logger import get_logger

parent = Path(__file__).parent
logger = get_logger(__name__)


def make_chem_pot_diag_from_mp(target: Union[Composition, str],
                               additional_elements: List[str] = None,
                               vertex_elements: List[str] = None,
                               atom_energy_yaml: Optional[str] = None):
    """Obtain the energies from Materials Project.

    When the atom_energy_yaml is provided, the total energies are aligned
    via atom energies.

    """
    properties = ["task_id", "full_formula", "final_energy"]
    target = target if isinstance(target, Composition) else Composition(target)
    elements = target.chemical_system.split("-")
    vertex_elements = vertex_elements or elements
    vertex_elements = [Element(e) for e in vertex_elements]
    if additional_elements:
        elements.extend(additional_elements)
    query = MpQuery(elements, properties=properties)
    comp_es = []
    if atom_energy_yaml:
        energies = loadfn(atom_energy_yaml)
        diff = {e: energies[e] - mp_energies[e] for e in elements}
    else:
        diff = None

    for m in query.materials:
        energy = m["final_energy"]
        if diff:
            for k, v in Composition(m["full_formula"]).as_dict().items():
                energy += diff[k] * v
        comp_es.append(CompositionEnergy(
            Composition(m["full_formula"]), energy, m["task_id"]))

    comp_es = remove_higher_energy_comp(comp_es)

    return ChemPotDiag(comp_es, target, vertex_elements)


def remove_higher_energy_comp(comp_energies: List[CompositionEnergy]):
    result = []
    for _, grouped_comp_energies in groupby(
            comp_energies, key=lambda x: x.composition.reduced_formula):
        result.append(min(list(grouped_comp_energies),
                          key=lambda y: y.abs_energy_per_atom))
    return result


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
