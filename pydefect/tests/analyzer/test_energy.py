# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.analyzer.energy import Energy
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def energy():
    return Energy(rel_energy=1.0, atom_io={"O": -1}, correction_energy=1.0,
                  is_shallow=False)


def test_energy_yaml(energy, tmpdir):
    expected_text = """rel_energy: 1.0
atom_io:
  O: -1
correction_energy: 1.0
is_shallow: False"""
    assert_yaml_roundtrip(energy, tmpdir, expected_text)


def test_energy_yaml2(tmpdir):
    energy = Energy(rel_energy=0.0, atom_io={})
    expected_text = """rel_energy: 1.0
atom_io:
"""
    assert_yaml_roundtrip(energy, tmpdir, expected_text)


