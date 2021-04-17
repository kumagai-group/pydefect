# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy

import pytest

from pydefect.analyzer.defect_energy import (
    DefectEnergy, CrossPoints, DefectEnergyInfo, DefectEnergies,
    DefectEnergySummary, ChargeEnergies)
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def defect_energy():
    return DefectEnergy(formation_energy=1.0,
                        energy_corrections={"1st order": 1.0, "alignment": 2.0},
                        is_shallow=False)


def test_energy_total_correction(defect_energy):
    assert defect_energy.total_correction == 3.0
    assert defect_energy.energy(with_correction=True) == 4.0


@pytest.fixture
def defect_energy_info(defect_energy):
    return DefectEnergyInfo(name="Va_O1", charge=1, atom_io={Element.O: -1},
                            defect_energy=defect_energy)


@pytest.fixture
def defect_energy_info2():
    energy = DefectEnergy(formation_energy=0.0,
                          energy_corrections={"no correction": 0.0})
    return DefectEnergyInfo(name="hole polaron", charge=1, atom_io={},
                            defect_energy=energy)


def test_energy_yaml(defect_energy_info, tmpdir):
    expected_text = """name: Va_O1
charge: 1
formation_energy: 1.0
atom_io:
  O: -1
energy_corrections:
  1st order: 1.0
  alignment: 2.0
is_shallow: False"""
    assert_yaml_roundtrip(defect_energy_info, tmpdir, expected_text)


def test_energy_yaml2(defect_energy_info2, tmpdir):
    expected_text = """name: hole polaron
charge: 1
formation_energy: 0.0
atom_io:
energy_corrections:
  no correction: 0.0
is_shallow: """
    assert_yaml_roundtrip(defect_energy_info2, tmpdir, expected_text)


@pytest.fixture
def defect_energies():
    return DefectEnergies(
        atom_io={Element.O: -1},
        charges=[0, 1, 2],
        defect_energies=[DefectEnergy(1.0, {"corr": 2.0}, is_shallow=False),
                         DefectEnergy(2.0, {"corr": 2.0}, is_shallow=False),
                         DefectEnergy(3.0, {"corr": 2.0}, is_shallow=True)])


@pytest.fixture
def defect_energy_summary(defect_energies):
    return DefectEnergySummary(
        title="MgAl2O4",
        defect_energies={"Va_O1": defect_energies},
        rel_chem_pots={"A": {Element.O: -1.0}, "B": {Element.O: -2.0}},
        cbm=2.0, supercell_vbm=-1.0, supercell_cbm=3.0)


def test_defect_energy_summary_charge_and_energies(
        defect_energy_summary, defect_energies):
    actual = defect_energy_summary.charge_and_energies(
        chem_pot_label="A", allow_shallow=True, with_correction=False)
    expected = ChargeEnergies(
        charge_energies=[(0, 0.0), (1, 1.0), (2, 2.0)], e_min=0.0, e_max=2.0)
    assert actual["Va_O"] == expected

    # change chem_pot_label
    actual = defect_energy_summary.charge_and_energies(
        chem_pot_label="B", allow_shallow=True, with_correction=False)
    expected = ChargeEnergies(
        charge_energies=[(0, -1.0), (1, 0.0), (2, 1.0)], e_min=0.0, e_max=2.0)
    assert actual["Va_O"] == expected

    # change allow_shallow
    actual = defect_energy_summary.charge_and_energies(
        chem_pot_label="A", allow_shallow=False, with_correction=False)
    expected = ChargeEnergies(
        charge_energies=[(0, 0.0), (1, 1.0)], e_min=0.0, e_max=2.0)
    assert actual["Va_O"] == expected

    # change with_correction
    actual = defect_energy_summary.charge_and_energies(
        chem_pot_label="A", allow_shallow=True, with_correction=True)
    expected = ChargeEnergies(
        charge_energies=[(0, 2.0), (1, 3.0), (2, 4.0)], e_min=0.0, e_max=2.0)
    assert actual["Va_O"] == expected

    # change e_min, e_max
    _defect_energy_summary = deepcopy(defect_energy_summary)
    _defect_energy_summary.e_min = -1.0
    _defect_energy_summary.e_max = 3.0
    actual = _defect_energy_summary.charge_and_energies(
        chem_pot_label="A", allow_shallow=True, with_correction=True)
    expected = ChargeEnergies(
        charge_energies=[(0, 2.0), (1, 3.0), (2, 4.0)], e_min=-1.0, e_max=3.0)
    assert actual["Va_O"] == expected


def test_defect_energy_summary_latexified_title(defect_energy_summary):
    actual = defect_energy_summary.latexified_title()
    expected = "MgAl$_{2}$O$_{4}$"
    assert actual == expected


@pytest.fixture
def charge_energies():
    return ChargeEnergies(charge_energies=[(0, 6.0), (1, 3.0)],
                          e_min=1.0, e_max=6.0)


def test_charge_energies_cross_points(charge_energies):
    actual = charge_energies.make_cross_point()
    expected = CrossPoints([[3.0, 6.0]], [[6.0, 6.0], [1.0, 4.0]])
    assert actual == expected


def test_charge_energies_transition_levels(charge_energies):
    actual = charge_energies.transition_levels
    expected = {(0, 1): 3.0}
    assert actual == expected


def test_charge_energies_pinning_level(charge_energies):
    assert charge_energies.pinning_level == (None, None)
    ce = ChargeEnergies(charge_energies=[(0, 6.0), (1, 3.0)],
                        e_min=-4.0, e_max=6.0)
    assert ce.pinning_level == ((-3.0, 1), None)


def test_charge_energies_energy_at_ef(charge_energies):
    assert charge_energies.energy_at_ef(ef=0.0) == (3.0, 1)
    assert charge_energies.energy_at_ef(ef=10.0) == (6.0, 0)


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
