# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CompositionEnergy
from pydefect.chem_pot_diag.make_chem_pot_diag import make_chem_pot_diag_from_mp
from pymatgen import Composition


@pytest.fixture
def mp_query(mocker):
    mock = mocker.patch("pydefect.chem_pot_diag.make_chem_pot_diag.MpQuery")
    mock.return_value.materials = \
        [{"full_formula": "O8", "final_energy": -39.58364375, "task_id": "mp-1"},
         {"full_formula": "Mg3", "final_energy": -4.79068775, "task_id": "mp-2"},
         {"full_formula": "Mg1O1", "final_energy": -11.96742144, "task_id": "mp-3"}]
    return mock


def test_make_chem_pot_diag_from_mp(mp_query):
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"))
    expected = ChemPotDiag({
        CompositionEnergy(Composition('O8'), -39.58364375, "mp-1"),
        CompositionEnergy(Composition('Mg3'), -4.79068775, "mp-2"),
        CompositionEnergy(Composition('Mg1O1'), -11.96742144, "mp-3")},
        target=Composition("MgO"))
    mp_query.assert_called_once_with(
        ["Mg", "O"], properties=["task_id", "full_formula", "final_energy"])
    assert actual == expected


def test_make_chem_pot_diag_from_mp_w_vise_functional(mp_query):
    to_datasets = Path(__file__).parent / "../../chem_pot_diag/datasets"
    pbesol = loadfn(str(to_datasets / "vise_pbesol_atom_energy.yaml"))
    mp = loadfn(str(to_datasets / "mp_atom_energy.yaml"))

    diff = {elem: pbesol[elem] - mp[elem] for elem in ["O", "Mg"]}
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"),
                                        vise_functional="pbesol")
    expected = ChemPotDiag({
        CompositionEnergy(Composition('O8'), -39.58364375 + diff["O"] * 8, "mp-1"),
        CompositionEnergy(Composition('Mg3'), -4.79068775 + diff["Mg"] * 3, "mp-2"),
        CompositionEnergy(Composition('Mg1O1'), -11.96742144 + diff["Mg"] + diff["O"], "mp-3")},
        target=Composition("MgO"))
    assert actual == expected


"""
TODO
- Make ChemPotDiag from MP data 
- Modify the MP total energies using the atom energy difference.
- Raise an error when the elements the MP does not support are given.

DONE
"""
