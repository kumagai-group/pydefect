# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from itertools import groupby
from pathlib import Path
from typing import List, Optional, Dict

from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergy, \
    CompositionEnergies
from pydefect.util.mp_tools import MpQuery
from pymatgen.core import Composition
from vise.atom_energies.atom_energy import mp_energies
from vise.util.logger import get_logger

parent = Path(__file__).parent
logger = get_logger(__name__)


def make_composition_energies_from_mp(elements: List[str],
                                      atom_energy_yaml: Optional[str] = None,
                                      ) -> CompositionEnergies:
    """Obtain the energies from Materials Project.

    When the atom_energy_yaml is provided, the total energies are aligned
    via atom energies.
    """
    properties = ["task_id", "full_formula", "final_energy"]
    query = MpQuery(elements, properties=properties)
    comp_es = {}
    if atom_energy_yaml:
        energies = loadfn(atom_energy_yaml)
        diff = {e: energies[e] - mp_energies[e] for e in elements}
    else:
        diff = {e: 0.0 for e in elements}

    for m in query.materials:
        key = Composition(m["full_formula"])
        energy = m["final_energy"]
        for k, v in key.as_dict().items():
            energy += diff[k] * v
        comp_es[key] = CompositionEnergy(energy, m["task_id"])
    comp_es = remove_higher_energy_comp(comp_es)
    return CompositionEnergies(comp_es)


def remove_higher_energy_comp(comp_energies: Dict[Composition, CompositionEnergy]):
    _l = [[k, v] for k, v in comp_energies.items()]
    result = {}
    for _, grouped_k_v in groupby(
            _l, key=lambda x: Composition(x[0]).reduced_formula):
        formula, comp_e = min(list(grouped_k_v),
                              key=lambda y: (y[1].energy
                                             / Composition(y[0]).num_atoms))
        result[formula] = comp_e
    return result
