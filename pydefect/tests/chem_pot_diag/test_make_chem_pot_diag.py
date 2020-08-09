# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag
from pydefect.chem_pot_diag.make_chem_pot_diag import make_chem_pot_diag_from_mp
from pymatgen import Composition


@pytest.fixture
def mp_query(mocker):
    mock = mocker.patch("pydefect.chem_pot_diag.make_chem_pot_diag.MpQuery")
    mock.return_value.materials = \
        [{"full_formula": "O8", "final_energy": -39.58364375},
         {"full_formula": "Mg3", "final_energy": -4.79068775},
         {"full_formula": "Mg1O1", "final_energy": -11.96742144}]
    return mock


def test_make_chem_pot_diag_from_mp(mp_query):
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"))
    expected = ChemPotDiag(energies={'O8': -39.58364375,
                                     'Mg3': -4.79068775,
                                     'Mg1O1': -11.96742144},
                           target=Composition("MgO"))
    mp_query.assert_called_once_with(
        ["Mg", "O"], properties=["task_id", "full_formula", "final_energy"])
    assert actual == expected


def test_make_chem_pot_diag_from_mp_w_vise_functional(mp_query):
    pbesol = loadfn("../../chem_pot_diag/datasets/vise_pbesol_atom_energy.yaml")
    mp = loadfn("../../chem_pot_diag/datasets/mp_atom_energy.yaml")

    diff = {elem: pbesol[elem] - mp[elem] for elem in ["O", "Mg"]}
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"),
                                        vise_functional="pbesol")
    expected = ChemPotDiag(
        energies={'O8': -39.58364375 + diff["O"] * 8,
                  'Mg3': -4.79068775 + diff["Mg"] * 3,
                  'Mg1O1': -11.96742144 + diff["Mg"] + diff["O"]},
        target=Composition("MgO"))
    assert actual == expected


"""
TODO
- Make ChemPotDiag from MP data 
- Modify the MP total energies using the atom energy difference.
- Raise an error when the elements the MP does not support are given.

DONE
"""
