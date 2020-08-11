# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path
from typing import List, Optional

from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CompositionEnergy
from pydefect.util.mp_tools import MpQuery
from pymatgen import Composition

pbesol_energies = loadfn(Path(__file__).parent / "datasets/vise_pbesol_atom_energy.yaml")
mp_energies = loadfn(Path(__file__).parent / "datasets/mp_atom_energy.yaml")


def make_chem_pot_diag_from_mp(elements: List[str],
                               target: Composition,
                               vise_functional: Optional[str] = None):
    """Obtain the energies from Materials Project."""
    properties = ["task_id", "full_formula", "final_energy"]
    query = MpQuery(elements, properties=properties)
    comp_es = set()
    for m in query.materials:
        energy = m["final_energy"]
        if vise_functional:
            diff = {elem: pbesol_energies[elem] - mp_energies[elem]
                    for elem in elements}
            for k, v in Composition(m["full_formula"]).as_dict().items():
                energy += diff[k] * v
        comp_es.add(CompositionEnergy(
            Composition(m["full_formula"]), energy, m["task_id"]))

    return ChemPotDiag(comp_es, target)


