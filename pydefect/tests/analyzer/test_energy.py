# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.analyzer.energy import EnergySummary, reservoir_energy, Energy
from pymatgen import Element
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def energy():
    return Energy(formation_energy=1.0,
                  atom_io={Element.O: -1},
                  correction_energy={"1st order": 1.0, "alignment-like": 2.0},
                  is_shallow=False)


def test_energy_total_correction(energy):
    assert energy.total_correction == 3.0


# def test_defect_formation_energy_wo_corr(energy, tmpdir):
#     actual = energy.formation_energy_wo_corr(abs_chem_pot={Element.O: 10.0})
#     assert actual == 1.0 - 10.0


# def test_defect_formation_energy_wo_correction(energy2, tmpdir):
#     energy2.correction_energy = {"a": 1.0}
#     actual = energy2.formation_energy_wo_corr(abs_chem_pot={})
#     assert actual == 0.0


@pytest.fixture
def energy_summary(energy):
    return EnergySummary(name="Va_O1", charge=1, energy=energy)


@pytest.fixture
def energy_summary2():
    energy = Energy(formation_energy=0.0, atom_io={},
                    correction_energy={"no correction": 0.0})
    return EnergySummary(name="hole polaron", charge=1, energy=energy)


def test_energy_yaml(energy_summary, tmpdir):
    expected_text = """name: Va_O1
charge: 1
formation_energy: 1.0
atom_io:
  O: -1
correction_energy:
  1st order: 1.0
  alignment-like: 2.0
is_shallow: False"""
    assert_yaml_roundtrip(energy_summary, tmpdir, expected_text)


def test_energy_yaml2(energy_summary2, tmpdir):
    expected_text = """name: hole polaron
charge: 1
formation_energy: 0.0
atom_io:
correction_energy:
  no correction: 0.0
is_shallow: """
    assert_yaml_roundtrip(energy_summary2, tmpdir, expected_text)


def test_reservoir_energy():
    num_atom_diff = {Element.He: 1, Element.Li: -1}
    abs_chem_pot = {Element.He: 10.0, Element.Li: 100.0, Element.H: 1000.0}
    actual = reservoir_energy(num_atom_diff, abs_chem_pot)
    assert actual == -10.0 + 100.0
