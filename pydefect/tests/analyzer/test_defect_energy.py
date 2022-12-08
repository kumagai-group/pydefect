# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.analyzer.defect_energy import (
    DefectEnergy, CrossPoints, DefectEnergyInfo, DefectEnergies,
    DefectEnergySummary, SingleChargeEnergies, ChargeEnergies)
from vise.tests.helpers.assertion import assert_yaml_roundtrip, \
    assert_json_roundtrip


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
    return DefectEnergyInfo(name="Va_O1", charge=1, atom_io={"O": -1},
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


def test_defect_energy_summary_json_roundtrip(defect_energy_summary, tmpdir):
    assert_json_roundtrip(defect_energy_summary, tmpdir)


def test_defect_energy_summary_str():
    defect_energies = DefectEnergies(
        atom_io={"O": -1},
        charges=[1, 2, 0],
        defect_energies=[
            DefectEnergy(2.0001, {"corr": 2.0001}, is_shallow=False),
            DefectEnergy(3.0001, {"corr": 2.0001}, is_shallow=True),
            DefectEnergy(1.0001, {"corr": 2.0001}, is_shallow=False),
        ])

    defect_energy_summary = DefectEnergySummary(
        title="MgAl2O4",
        defect_energies={"Va_O1": defect_energies},
        rel_chem_pots={"A": {"O": -1.000000001}, "B": {"O": -2.000000001}},
        cbm=2.000000001, supercell_vbm=-1.000000001, supercell_cbm=3.000000001)

    actual = defect_energy_summary.__str__()
    expected = """title: MgAl₂O₄
rel_chem_pots:
 -A O: -1.00
 -B O: -2.00
vbm: 0.00, cbm: 2.00, supercell vbm: -1.00, supercell_cbm: 3.00

name    atom_io      charge    energy    correction  is_shallow
------  ---------  --------  --------  ------------  ------------
Va_O1   O: -1             0     1.000         2.000  False
                          1     2.000         2.000  False
                          2     3.000         2.000  True"""
    assert actual == expected


def test_defect_energy_summary_charge_energies(defect_energy_summary):
    actual = defect_energy_summary.charge_energies(
        chem_pot_label="A", allow_shallow=True, with_corrections=False, e_range=(0.0, 2.0))
    expected = ChargeEnergies(
        charge_energies_dict={"Va_O": SingleChargeEnergies([(0, 0.0), (1, 1.0), (2, 2.0)])},
        e_min=0.0, e_max=2.0)
    assert actual == expected

    # change chem_pot_label
    actual = defect_energy_summary.charge_energies(
        chem_pot_label="B", allow_shallow=True, with_corrections=False, e_range=(0.0, 2.0))
    expected = ChargeEnergies(
        charge_energies_dict={"Va_O": SingleChargeEnergies([(0, -1.0), (1, 0.0), (2, 1.0)])},
        e_min=0.0, e_max=2.0)
    assert actual == expected

    # change allow_shallow
    actual = defect_energy_summary.charge_energies(
        chem_pot_label="A", allow_shallow=False, with_corrections=False, e_range=(0.0, 2.0))
    expected = ChargeEnergies(
        charge_energies_dict={"Va_O": SingleChargeEnergies([(0, 0.0), (1, 1.0)])},
        e_min=0.0, e_max=2.0)
    assert actual == expected

    # change with_correction
    actual = defect_energy_summary.charge_energies(
        chem_pot_label="A", allow_shallow=True, with_corrections=True, e_range=(0.0, 2.0))
    expected = ChargeEnergies(
        charge_energies_dict={"Va_O": SingleChargeEnergies([(0, 2.0), (1, 3.0), (2, 4.0)])},
        e_min=0.0, e_max=2.0)
    assert actual == expected


def test_cross_points_e_max_energies_dicts():
    charge_energies = ChargeEnergies(
        charge_energies_dict={"Va_Mg": SingleChargeEnergies([(0, 2.0)]),
                              "Va_O": SingleChargeEnergies([(1, 1.0)])},
        e_min=0.0, e_max=2.0)
    actual = charge_energies.cross_point_dicts
    expected = {"Va_Mg": CrossPoints(inner_cross_points=[],
                                     boundary_points=[[0.0, 2.0], [2.0, 2.0]]),
                "Va_O": CrossPoints(inner_cross_points=[],
                                    boundary_points=[[0.0, 1.0], [2.0, 3.0]])}
    assert actual == expected

    actual = charge_energies.e_min_max_energies_dict
    expected = {"Va_Mg": [[2.0, 2.0]], "Va_O": [[1.0, 3.0]]}
    assert actual == expected


def test_defect_energy_summary_latexified_title(defect_energy_summary):
    actual = defect_energy_summary.latexified_title
    expected = "MgAl$_{2}$O$_{4}$"
    assert actual == expected


@pytest.fixture
def single_charge_energies():
    return SingleChargeEnergies([(0, 6.0), (1, 3.0)])


@pytest.fixture
def charge_energies(single_charge_energies):
    return ChargeEnergies({"Va_O": single_charge_energies},
                          e_min=1.0, e_max=6.0)


def test_charge_energies_cross_points(charge_energies):
    actual = charge_energies.cross_point_dicts
    expected = {"Va_O": CrossPoints([[3.0, 6.0]], [[6.0, 6.0], [1.0, 4.0]])}
    assert actual == expected


def test_charge_energies_energy_range(charge_energies):
    actual = charge_energies.energy_range(0.2)
    expected = [3.8, 6.2]
    assert actual == expected


def test_charge_energies_pinning_level(single_charge_energies):
    assert single_charge_energies.pinning_level(1.0, 6.0) == (None, None)

    ce_1 = SingleChargeEnergies(charge_energies=[(0, 6.0), (1, 3.0)])
    assert ce_1.pinning_level(-4.0, 6.0) == ((-3.0, 1), None)

    ce_2 = SingleChargeEnergies(charge_energies=[(1, -1.0), (2, 0.0)])
    assert ce_2.pinning_level(0, 2.0) == ((1.0, 1), None)


def test_charge_energies_energy_at_ef(single_charge_energies):
    assert single_charge_energies.energy_at_ef(ef=0.0) == (3.0, 1)
    assert single_charge_energies.energy_at_ef(ef=10.0) == (6.0, 0)


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


