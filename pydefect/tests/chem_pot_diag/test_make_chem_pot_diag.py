# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CompositionEnergy
from pydefect.chem_pot_diag.make_chem_pot_diag import \
    make_chem_pot_diag_from_mp, remove_higher_energy_comp
from pymatgen import Composition


@pytest.fixture
def mp_query(mocker):
    mock = mocker.patch("pydefect.chem_pot_diag.make_chem_pot_diag.MpQuery")
    mock.return_value.materials = \
        [{"full_formula": "O8", "final_energy": -39.58364375, "task_id": "mp-1"},
         {"full_formula": "Mg3", "final_energy": -4.79068775, "task_id": "mp-2"},
         {"full_formula": "Mg1O1", "final_energy": -11.96742144, "task_id": "mp-3"}]
    return mock


@pytest.fixture
def cpd():
    return ChemPotDiag({
        CompositionEnergy(Composition('O8'), -39.58364375, "mp-1"),
        CompositionEnergy(Composition('Mg3'), -4.79068775, "mp-2"),
        CompositionEnergy(Composition('Mg1O1'), -11.96742144, "mp-3")},
        target=Composition("MgO"))


to_datasets = Path(__file__).parent / "../../chem_pot_diag/datasets"
pbesol = loadfn(str(to_datasets / "vise_pbesol_atom_energy.yaml"))
mp = loadfn(str(to_datasets / "mp_atom_energy.yaml"))

diff = {elem: pbesol[elem] - mp[elem] for elem in ["O", "Mg"]}


@pytest.fixture
def cpd_corr():
    return ChemPotDiag({
        CompositionEnergy(Composition('O8'), -39.58364375 + diff["O"] * 8, "mp-1"),
        CompositionEnergy(Composition('Mg3'), -4.79068775 + diff["Mg"] * 3, "mp-2"),
        CompositionEnergy(Composition('Mg1O1'), -11.96742144 + diff["Mg"] + diff["O"], "mp-3")},
        target=Composition("MgO"))


def test_make_chem_pot_diag_from_mp(mp_query, cpd):
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"))
    mp_query.assert_called_once_with(
        ["Mg", "O"], properties=["task_id", "full_formula", "final_energy"])
    assert actual == cpd


def test_make_chem_pot_diag_from_mp_yaml(mp_query, tmpdir, cpd_corr):
    tmpdir.chdir()
    tmpdir.join("tmp.yaml").write("""O:  -1.61154565
Mg: -0.00912097""")
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"),
                                        atom_energy_yaml="tmp.yaml")
    assert actual == cpd_corr


def test_make_chem_pot_diag_from_mp_w_vise_functional(mp_query, cpd_corr):
    actual = make_chem_pot_diag_from_mp(elements=["Mg", "O"],
                                        target=Composition("MgO"),
                                        atom_energy_yaml="pbesol")
    assert actual == cpd_corr


def test_remove_higher_energy_comp():
    comp_es = {CompositionEnergy(Composition('O2'), -2.1, "mp-1"),
               CompositionEnergy(Composition('O4'), -4.0, "mp-2"),
               CompositionEnergy(Composition('Mg1'), -10, "mp-3")}
    actual = remove_higher_energy_comp(comp_es)
    expected = {CompositionEnergy(Composition('O2'), -2.1, "mp-1"),
                CompositionEnergy(Composition('Mg1'), -10, "mp-3")}
    assert actual == expected

"""
TODO
- Raise an error when the elements the MP does not support are given.

DONE
"""
