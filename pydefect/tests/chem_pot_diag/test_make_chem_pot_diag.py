# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CompositionEnergy
from pydefect.chem_pot_diag.make_chem_pot_diag import \
    make_chem_pot_diag_from_mp, remove_higher_energy_comp
from pymatgen import Composition, Element


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
    return ChemPotDiag([
        CompositionEnergy(Composition('O8'), -39.58364375, "mp-1"),
        CompositionEnergy(Composition('Mg3'), -4.79068775, "mp-2"),
        CompositionEnergy(Composition('Mg1O1'), -11.96742144, "mp-3")],
        target=Composition("MgO"))
# mp
# O:  -1.54794890
# Mg: -0.01731593

diff = {"Mg": -0.00912097 - -0.01731593, "O": -1.61154565 - -1.54794890}


@pytest.fixture
def cpd_corr():
    return ChemPotDiag([
        CompositionEnergy(Composition('O8'), -39.58364375 + diff["O"] * 8, "mp-1"),
        CompositionEnergy(Composition('Mg3'), -4.79068775 + diff["Mg"] * 3, "mp-2"),
        CompositionEnergy(Composition('Mg1O1'), -11.96742144 + diff["Mg"] + diff["O"], "mp-3")],
        target=Composition("MgO"))


@pytest.mark.parametrize("composition", [Composition("MgO"), "MgO"])
def test_make_chem_pot_diag_from_mp(composition, mp_query, cpd):
    actual = make_chem_pot_diag_from_mp(target=composition)
    mp_query.assert_called_once_with(
        ["Mg", "O"], properties=["task_id", "full_formula", "final_energy"])
    assert actual == cpd
    assert actual.target == cpd.target
    assert actual.vertex_elements == cpd.vertex_elements


def test_make_chem_pot_diag_from_mp_yaml(mp_query, tmpdir, cpd_corr):
    tmpdir.chdir()
    tmpdir.join("tmp.yaml").write("""O:  -1.61154565
Mg: -0.00912097""")
    actual = make_chem_pot_diag_from_mp(target=Composition("MgO"),
                                        atom_energy_yaml="tmp.yaml")
    assert actual == cpd_corr


def test_make_chem_pot_diag_from_mp_additional_elem(mp_query, tmpdir, cpd_corr):
    tmpdir.chdir()
    tmpdir.join("tmp.yaml").write("""O:  -1.61154565
Mg: -0.00912097""")
    actual = make_chem_pot_diag_from_mp(target=Composition("Mg"),
                                        additional_elements=["O"],
                                        atom_energy_yaml="tmp.yaml",
                                        vertex_elements=["Mg", "O"])
    assert actual.vertex_elements == [Element.Mg, Element.O]


def test_remove_higher_energy_comp():
    comp_es = [CompositionEnergy(Composition('O2'), -2.1, "mp-1"),
               CompositionEnergy(Composition('O4'), -4.0, "mp-2"),
               CompositionEnergy(Composition('Mg1'), -10, "mp-3")]
    actual = remove_higher_energy_comp(comp_es)
    expected = [CompositionEnergy(Composition('O2'), -2.1, "mp-1"),
                CompositionEnergy(Composition('Mg1'), -10, "mp-3")]
    assert actual == expected

