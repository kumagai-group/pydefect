# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path
from typing import List, Optional

import fire
from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag
from pydefect.util.mp_tools import MpQuery
from pymatgen import Composition


def make_chem_pot_diag_from_mp(elements: List[str],
                               target: Composition,
                               vise_functional: Optional[str] = None):
    """Obtain the energies from Materials Project."""
    properties = ["task_id", "full_formula", "final_energy"]
    query = MpQuery(elements, properties=properties)
    energies = {m["full_formula"]: m["final_energy"] for m in query.materials}

    if vise_functional:
        pbesol_energies = loadfn(
            Path(__file__).parent / "datasets/vise_pbesol_atom_energy.yaml")
        mp_energies = loadfn(
            Path(__file__).parent / "datasets/mp_atom_energy.yaml")
        diff = {elem: pbesol_energies[elem] - mp_energies[elem]
                for elem in elements}
        for formula, energy in energies.items():
            for k, v in Composition(formula).as_dict().items():
                energies[formula] += diff[k] * v

    cpd = ChemPotDiag(energies, target)
    cpd.to_json_file("chem_pot_diag_mp.json")
    return cpd


if __name__ == '__main__':
    fire.Fire(make_chem_pot_diag_from_mp)

