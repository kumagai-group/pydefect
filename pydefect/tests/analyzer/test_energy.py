# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.analyzer.energy import Energy, reservoir_energy
from pymatgen import Element
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def energy():
    return Energy(name="Va_O1", charge=1, rel_energy=1.0, atom_io={Element.O: -1},
                  correction_energy={"1st order": 1.0, "alignment-like": 2.0},
                  is_shallow=False)


@pytest.fixture
def energy2():
    return Energy(name="hole polaron", charge=1, rel_energy=0.0, atom_io={})


def test_energy_yaml(energy, tmpdir):
    expected_text = """name: Va_O1
charge: 1
rel_energy: 1.0
atom_io:
  O: -1
correction_energy: 
  1st order: 1.0
  alignment-like: 2.0
is_shallow: False"""
    assert_yaml_roundtrip(energy, tmpdir, expected_text)


def test_energy_yaml2(energy2, tmpdir):
    expected_text = """"name: Va_O1
charge: 1
rel_energy: 1.0
atom_io:
"""
    assert_yaml_roundtrip(energy2, tmpdir, expected_text)


def test_energy_total_correction(energy, energy2, tmpdir):
    assert energy.total_correction == 3.0
    assert energy2.total_correction is None


def test_defect_formation_energy(energy, tmpdir):
    actual = energy.defect_formation_energy(abs_chem_pot={Element.O: 10.0})
    assert actual == 1.0 - 10.0 + 3.0


def test_defect_formation_energy_wo_correction(energy2, tmpdir):
    with pytest.raises(ValueError):
        energy2.defect_formation_energy(abs_chem_pot={})
    energy2.correction_energy = {"a": 1.0}
    actual = energy2.defect_formation_energy(abs_chem_pot={})
    assert actual == 0.0 + 1.0


def test_reservoir_energy():
    num_atom_diff = {Element.He: 1, Element.Li: -1}
    abs_chem_pot = {Element.He: 10.0, Element.Li: 100.0, Element.H: 1000.0}
    actual = reservoir_energy(num_atom_diff, abs_chem_pot)
    assert actual == -10.0 + 100.0
