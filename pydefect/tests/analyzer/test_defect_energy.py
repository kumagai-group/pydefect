# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.analyzer.defect_energy import (
    DefectEnergy, CrossPoints, defect_mpl_name,
    slide_energy, sanitize_defect_energies_for_plot,
    defect_plotly_name, make_defect_energies, DefectEnergyInfo,
    reservoir_energy)
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def defect_energy():
    return DefectEnergy(formation_energy=1.0,
                        correction_energy={"1st order": 1.0, "alignment": 2.0},
                        is_shallow=False)


def test_energy_total_correction(defect_energy):
    assert defect_energy.total_correction == 3.0


@pytest.fixture
def defect_energy_info(defect_energy):
    return DefectEnergyInfo(name="Va_O1", charge=1, atom_io={Element.O: -1},
                            energy=defect_energy)


@pytest.fixture
def defect_energy_info2():
    energy = DefectEnergy(formation_energy=0.0,
                          correction_energy={"no correction": 0.0})
    return DefectEnergyInfo(name="hole polaron", charge=1, atom_io={}, energy=energy)


def test_energy_yaml(defect_energy_info, tmpdir):
    expected_text = """name: Va_O1
charge: 1
formation_energy: 1.0
atom_io:
  O: -1
correction_energy:
  1st order: 1.0
  alignment: 2.0
is_shallow: False"""
    assert_yaml_roundtrip(defect_energy_info, tmpdir, expected_text)


def test_energy_yaml2(defect_energy_info2, tmpdir):
    expected_text = """name: hole polaron
charge: 1
formation_energy: 0.0
atom_io:
correction_energy:
  no correction: 0.0
is_shallow: """
    assert_yaml_roundtrip(defect_energy_info2, tmpdir, expected_text)


def test_reservoir_energy():
    num_atom_diff = {Element.He: 1, Element.Li: -1}
    abs_chem_pot = {Element.He: 10.0, Element.Li: 100.0, Element.H: 1000.0}
    actual = reservoir_energy(num_atom_diff, abs_chem_pot)
    assert actual == -10.0 + 100.0


@pytest.fixture
def defect_energies():
    return DefectEnergy(name="Va_O1",
                        charges=[0, 1, 2],
                        energies=[4, 2, -4],
                        corrections=[2, 1, 0])


def test_defect_energies_cross_points(defect_energies):
    actual = defect_energies.cross_points(1, 6)
    expected = CrossPoints([[5.0, 6.0]], [[1.0, -2.0], [6.0, 6.0]])
    assert actual == expected

    actual = defect_energies.cross_points(1, 6, 1)
    expected = CrossPoints([[4.0, 6.0]], [[0.0, -2.0], [5.0, 6.0]])
    assert actual == expected


def test_defect_energies_transition_levels(defect_energies):
    actual = defect_energies.transition_levels(base_e=0.0)
    expected = {(0, 1): 3.0, (0, 2): 5.0, (1, 2): 7.0}
    assert actual == expected

    actual = defect_energies.transition_levels(base_e=0.1)
    expected = {(0, 1): 2.9, (0, 2): 4.9, (1, 2): 6.9}
    assert actual == expected


def test_defect_energies_str(defect_energies):
    expected = """     Va_O1    0       4.0000       2.0000
     Va_O1    1       2.0000       1.0000
     Va_O1    2      -4.0000       0.0000"""
    assert str(defect_energies) == expected


def test_stable_charges(defect_energies):
    actual = defect_energies.stable_charges(ef_min=4.9, ef_max=5.1)
    assert actual == {0, 2}
    actual = defect_energies.stable_charges(ef_min=0.0, ef_max=0.1)
    assert actual == {2}


def test_pinning(defect_energies):
    assert defect_energies.pinning_level() == ((2.0, 2), (float("inf"), None))
    assert defect_energies.pinning_level(base_e=1.0) == ((1.0, 2), (float("inf"), None))


def test_energy_at_ef(defect_energies):
    assert defect_energies.energy_at_ef(ef=0.0) == (-4.0, 2)
    assert defect_energies.energy_at_ef(ef=10.0) == (6.0, 0)


@pytest.fixture
def cross_points():
    inner_cross_points = [[2, 30], [3, 40]]
    boundary_points = [[1, 10], [4, 40]]
    return CrossPoints(inner_cross_points, boundary_points)


def test_cross_points(cross_points):
    assert cross_points.all_sorted_points == [[1, 10], [2, 30], [3, 40], [4, 40]]
    assert cross_points.t_all_sorted_points == [[1, 2, 3, 4], [10, 30, 40, 40]]
    assert cross_points.t_inner_cross_points == [[2, 3], [30, 40]]
    assert cross_points.t_boundary_points == [[1, 4], [10, 40]]
    assert cross_points.charges == [20, 10, 0]
    assert cross_points.charge_list == [(None, 20), (20, 10), (10, 0), (0, None)]
    assert cross_points.annotated_charge_positions == \
           {20: [1.5, 20], 10: [2.5, 35], 0: [3.5, 40]}


def test_cross_points_str(cross_points):
    expected = """      1.0000      10.0000
      2.0000      30.0000
      3.0000      40.0000
      4.0000      40.0000"""
    assert str(cross_points) == expected


def test_make_defect_energies():
    defect_energies = [
        DefectEnergyInfo(name="Va_Mg1", charge=0, rel_energy=-2.0, atom_io={Element.Mg: -1},
                         correction_energy={"PC correction": 2.0}, is_shallow=False),
        DefectEnergyInfo(name="Va_Mg1", charge=-1, rel_energy=-14.0, atom_io={Element.Mg: -1},
                         correction_energy={"PC correction": 4.0}, is_shallow=False),
        DefectEnergyInfo(name="Va_Mg1", charge=-2, rel_energy=-14.0, atom_io={Element.Mg: -1},
                         correction_energy={"PC correction": 6.0}, is_shallow=True),

        DefectEnergyInfo(name="Va_O1", charge=0, rel_energy=-1.0, atom_io={Element.O: -1},
                         correction_energy={"PC correction": 1.0}, is_shallow=False),
        DefectEnergyInfo(name="Va_O1", charge=1, rel_energy=7.0, atom_io={Element.O: -1},
                         correction_energy={"PC correction": 3.0}, is_shallow=False),
        DefectEnergyInfo(name="Va_O1", charge=2, rel_energy=15.0, atom_io={Element.O: -1},
                         correction_energy={"PC correction": 5.0}, is_shallow=False),
    ]
    abs_chem_pot = {Element.Mg: 5.0, Element.O: 3.0}
    actual = make_defect_energies(defect_energies, abs_chem_pot, allow_shallow=True)
    expected = [DefectEnergy("Va_Mg1", [0, -1, -2], [-7.0, -19.0, -19.0], [2.0, 4.0, 6.0]),
                DefectEnergy("Va_O1", [0, 1, 2], [-4.0, 4.0, 12.0], [1.0, 3.0, 5.0])]
    assert actual == expected

    actual = make_defect_energies(defect_energies, abs_chem_pot, allow_shallow=False)
    expected = [DefectEnergy("Va_Mg1", [0, -1], [-7.0, -19.0], [2.0, 4.0]),
                DefectEnergy("Va_O1", [0, 1, 2], [-4.0, 4.0, 12.0], [1.0, 3.0, 5.0])]
    assert actual == expected


def test_defect_mpl_name():
    assert defect_mpl_name(name="Va_O1") == "$V_{{\\rm O}1}$"
    assert defect_mpl_name(name="Mg_i1") == "${\\rm Mg}_{i1}$"


def test_defect_plotly_name():
    assert defect_plotly_name(name="Va_O1") == "<i>V</i><sub>O1</sub>"


def test_generate_defect_energies():
    energies = [DefectEnergy("Va_Mg1", [0], [0], [0]),
                DefectEnergy("Va_O1", [0], [0], [0]),
                DefectEnergy("Va_O2", [0], [0], [0]),
                DefectEnergy("Mg_i1", [0], [0], [0]),
                DefectEnergy("O_i1", [0], [0], [0])]
    actual = sanitize_defect_energies_for_plot(energies)
    expected = [DefectEnergy("$V_{{\\rm Mg}}$", [0], [0], [0]),
                DefectEnergy("$V_{{\\rm O}1}$", [0], [0], [0]),
                DefectEnergy("$V_{{\\rm O}2}$", [0], [0], [0]),
                DefectEnergy("${\\rm Mg}_{i}$", [0], [0], [0]),
                DefectEnergy("${\\rm O}_{i}$", [0], [0], [0])]
    assert actual == expected

    energies = [DefectEnergy("Mg_i1", [0], [0], [0]),
                DefectEnergy("Mg_i2", [0], [0], [0]),]
    actual = sanitize_defect_energies_for_plot(energies)
    expected = [DefectEnergy("${\\rm Mg}_{i1}$", [0], [0], [0]),
                DefectEnergy("${\\rm Mg}_{i2}$", [0], [0], [0])]
    assert actual == expected

    energies = [DefectEnergy("Mg_i1", [0], [0], [0]),
                DefectEnergy("Mg_i2", [0], [0], [0]),]
    actual = sanitize_defect_energies_for_plot(energies, for_plotly=True)
    expected = [DefectEnergy("Mg<sub>i1</sub>", [0], [0], [0]),
                DefectEnergy("Mg<sub>i2</sub>", [0], [0], [0])]
    assert actual == expected


def test_slide_energy():
    energies = [DefectEnergy("Va_Mg1", [0], [0], [0]),
                DefectEnergy("Va_O1", [1], [0], [0]),
                DefectEnergy("Va_O2", [2], [0], [0]),
                DefectEnergy("Mg_i1", [0], [0], [0])]
    actual = slide_energy(energies, 1.0)
    expected = [DefectEnergy("Va_Mg1", [0], [0], [0]),
                DefectEnergy("Va_O1", [1], [1.0], [0]),
                DefectEnergy("Va_O2", [2], [2.0], [0]),
                DefectEnergy("Mg_i1", [0], [0], [0])]
    assert actual == expected

    """
    TODO
    - Evaluate the crossing points at given Fermi level range.
    """