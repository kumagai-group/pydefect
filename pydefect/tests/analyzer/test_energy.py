# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.analyzer.energy import Energy
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
