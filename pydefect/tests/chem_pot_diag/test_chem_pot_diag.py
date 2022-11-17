# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy

import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, \
    CompositionEnergy, CompositionEnergies, StandardEnergies, \
    RelativeEnergies, ChemPotDiagMaker, TargetVertices, TargetVertex, \
    target_element_chem_pot, change_element_sequence, UnstableTargetError
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.core import Composition
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def composition_energies():
    return CompositionEnergies(
        {Composition("H"): CompositionEnergy(0.0, "a"),
         Composition("H4O2"): CompositionEnergy(-4.0, "c"),
         Composition("O"): CompositionEnergy(1.0, "b"),
         Composition("Cl"): CompositionEnergy(12.0, "f"),
         Composition("O2Cl"): CompositionEnergy(3.0, "e"),
         Composition("O2Cl2"): CompositionEnergy(6.0, "d")})


def test_composition_energies_yaml_roundtrip(composition_energies, tmpdir):
    expected_text = """Cl1:
  energy: 12.0
  source: f
H1:
  energy: 0.0
  source: a
H4O2:
  energy: -4.0
  source: c
O1:
  energy: 1.0
  source: b
O2Cl1:
  energy: 3.0
  source: e
O2Cl2:
  energy: 6.0
  source: d
"""
    assert_yaml_roundtrip(composition_energies, tmpdir, expected_text,
                          compare_dict=False, compare_items=False)


def test_composition_energies_from_dict():
    actual = CompositionEnergies.from_dict({"Mg": -1.0})
    expected = CompositionEnergies({Composition("Mg"): CompositionEnergy(-1.0)})
    assert actual == expected
    assert actual.std_rel_energies == ({'Mg': -1.0}, {"Mg": 0.0})


def test_composition_energies_elements(composition_energies):
    assert composition_energies.elements == ["Cl", "H", "O"]


def test_std_rel_energies(composition_energies):
    actual = composition_energies.std_rel_energies
    expected_ref = StandardEnergies({"H": 0.0, "O": 1.0, "Cl": 12.0})
    expected_rel = RelativeEnergies({"H": 0.0, "O": 0.0, "Cl": 0.0,
                                     'ClO': -5.0, 'ClO2': -3.666666666666667,
                                     'H2O': -1.0})
    assert actual[0] == expected_ref
    assert actual[1] == expected_rel


@pytest.fixture
def standard_energies():
    return StandardEnergies({"Mg": 1.0, "O": 2.0})


def test_std_energies(standard_energies, tmpdir):
    expected_text = """Mg: 1.0
O: 2.0
"""
    assert_yaml_roundtrip(standard_energies, tmpdir, expected_text,
                          compare_dict=False, compare_items=False)


def test_target_element_chem_pot():
    actual = target_element_chem_pot("MgO3",
                                     energy_per_atom=-10.0,
                                     target_element="Mg",
                                     other_elem_chem_pot={"O": -4.0})
    assert actual == (-10.0 - (-4.0 * 3 / 4)) * 4


@pytest.fixture
def relative_energies():
    return RelativeEnergies({"Mg": 0.0, "MgO2": -3.0, "MgAlO2": -1.0,
                             "Al": 0.0, "O": 0.0})


def test_relative_energies_yaml(relative_energies, tmpdir):
    expected_text = """Al: 0.0
Mg: 0.0
MgAlO2: -1.0
MgO2: -3.0
O: 0.0
"""
    assert_yaml_roundtrip(relative_energies, tmpdir, expected_text,
                          compare_dict=False, compare_items=False)


def test_relative_energies_element_set(relative_energies):
    assert relative_energies.all_element_set == {"Mg", "Al", "O"}


def test_relative_energies_related_comp_energies(relative_energies):
    actual = relative_energies.host_composition_energies(elements=["Mg", "O"])
    expected = {"Mg": 0.0, "MgO2": -3.0, "O": 0.0}
    assert actual == expected


def test_relative_energies_comp_energies_w_element(relative_energies):
    actual = relative_energies.comp_energies_with_element(element="Al")
    expected = {"Al": 0.0, "MgAlO2": -1.0}
    assert actual == expected


def test_relative_energies_phase_diagram():
    rel_energies = RelativeEnergies({"Mg": 0.0, "MgO": -20.0, "MgO3": -4.0,
                                     "O": 0.0})
    expected = {PDEntry(Composition("MgO3"), energy=-4.0):
                    ({PDEntry(Composition("MgO"), -20.0): 0.5,
                      PDEntry(Composition("O"), 0.0): 0.5},
                     4.0)}
    assert rel_energies.unstable_compounds == expected
    expected = """|        |   composition | E above hull   | decompose to (ratio)   |
|--------+---------------+----------------+------------------------|
| Mg1 O3 |             4 | Mg1 O1 (0.500) | O1 (0.500)             |"""
    assert rel_energies.unstable_comp_info == expected


@pytest.fixture
def cpd_maker():
    return ChemPotDiagMaker(RelativeEnergies({"MgO2": -3.0, "MgAlO2": -2.5}),
                            elements=["Mg", "O"])


def test_cpd_maker_chem_pot_diag(cpd_maker):
    actual = cpd_maker.chem_pot_diag
    min_val = -9.0 * 1.1
    expected = ChemPotDiag(vertex_elements=["Mg", "O"],
                           polygons={'Mg': [[0.0, min_val], [0.0, -4.5]],
                                     'MgO2': [[0.0, -4.5], [-9.0, 0.0]],
                                     'O': [[min_val, 0.0], [-9.0, 0.0]]})
    assert actual == expected


def test_chem_pot_diag_maker_for_unstable_target_compound():
    cpd_maker = ChemPotDiagMaker(RelativeEnergies({"MgO": 0.1}),
                                 elements=["Mg", "O"],
                                 target="MgO")
    with pytest.raises(UnstableTargetError):
        cpd_maker.chem_pot_diag


@pytest.fixture
def cpd():
    min_val = -9.0 * 1.1
    return ChemPotDiag(
        vertex_elements=["Mg", "O"],
        polygons={'Mg': [[0.0, min_val], [0.0, -4.5]],
                  'MgO2': [[0.0, -4.5], [-9.0, 0.0]],
                  'O': [[min_val, 0.0], [-9.0, 0.0]]},
        target="MgO2",
        target_vertices_dict={
            "A": TargetVertex({"Mg": 0.0, "O": -4.5, "Al": -1.0},
                              ["Mg"], ["MgAlO2"]),
            "B": TargetVertex({"Mg": -9.0, "O": 0.0, "Al": -1.0},
                              ["O"], ["MgAlO2"])})


@pytest.fixture
def target_vertices():
    vertex_a = TargetVertex(chem_pot={"Mg": 0.0, "O": -4.5, "Al": -1.0},
                            competing_phases=["Mg"],
                            impurity_phases=["MgAlO2"])
    vertex_b = TargetVertex(chem_pot={"Mg": -9.0, "O": 0.0, "Al": -1.0},
                            competing_phases=["O"],
                            impurity_phases=["MgAlO2"])
    return TargetVertices(target="MgO2",
                          vertices={"A": vertex_a, "B": vertex_b})


def test_cpd_maker_chem_pot_diag_with_target(cpd_maker, cpd, target_vertices):
    cpd_maker_2 = copy(cpd_maker)
    cpd_maker_2.target = "MgO2"
    actual = cpd_maker_2.chem_pot_diag
    assert actual == cpd


def test_target_vertices(target_vertices, tmpdir):
    tmpdir.chdir()
    expected_text = """target: MgO2
A:
  chem_pot:
    Al: -1.0
    Mg: 0.0
    O: -4.5
  competing_phases:
  - Mg
  impurity_phases:
  - MgAlO2
B:
  chem_pot:
    Al: -1.0
    Mg: -9.0
    O: 0.0
  competing_phases:
  - O
  impurity_phases:
  - MgAlO2
"""
    assert_yaml_roundtrip(target_vertices, tmpdir, expected_text)


def test_cpd(cpd):
    assert cpd.target_coords == {"A": [0.0, -4.5], "B": [-9.0, 00]}
    assert cpd.dim == 2
    assert cpd.comp_centers == {"Mg": [0.0, (-9.9 - 4.5) / 2],
                                "MgO2": [-4.5, -2.25],
                                "O": [(-9.9 - 9.0) / 2, 0.0]}
    assert cpd.min_range == -9.9
    assert cpd.chemical_system == "Mg-O"


def test_cpd_change_element_sequence(cpd):
    min_val = -9.0 * 1.1
    actual = change_element_sequence(cpd, element_sequence=["O", "Mg"])
    assert actual.vertex_elements == ["O", "Mg"]
    assert actual.polygons == {'Mg': [[min_val, 0.0], [-4.5, 0.0]],
                               'MgO2': [[-4.5, 0.0], [0.0, -9.0]],
                               'O': [[0.0, min_val], [0.0, -9.0]]}

    actual = change_element_sequence(actual)
    assert actual == cpd

    actual.target = None
    assert change_element_sequence(actual).polygons == cpd.polygons
