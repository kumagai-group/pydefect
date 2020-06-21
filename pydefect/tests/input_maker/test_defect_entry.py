# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import numpy as np
import pytest
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.tests.helpers.assertion import assert_json_roundtrip, \
    assert_msonable
from pymatgen import IStructure, Lattice

# "H" at [0.0, 0.0, 0.0] is removed here.
coords = \
    [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
     [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.5]]
perturbed_coords = \
    [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
     [0.0, 0.0, 0.501], [0.0, 0.501, 0.0], [0.501, 0.0, 0.0], [0.5, 0.5, 0.5]]

rocksalt = IStructure(
    Lattice.cubic(10.0), ["H"] * 3 + ["He"] * 4, coords)
perturbed_rocksalt = IStructure(
    Lattice.cubic(10.0), ["H"] * 3 + ["He"] * 4, perturbed_coords)


@pytest.fixture
def defect_entry():
    return DefectEntry(name="Va_O1",
                       charge=1,
                       structure=rocksalt,
                       perturbed_structure=perturbed_rocksalt,
                       site_symmetry="m-3m",
                       defect_center=(0, 0, 0))


def test_hashable(defect_entry):
    d = {defect_entry: 1}


def test_perturbed_site_indices(defect_entry):
    assert defect_entry.perturbed_site_indices == [3, 4, 5]


def test_full_name(defect_entry):
    assert defect_entry.full_name == "Va_O1_1"


def test_anchor_atom_index(defect_entry):
    assert defect_entry.anchor_atom_index == 6


def test_anchor_atom_coords(defect_entry):
    np.testing.assert_array_almost_equal(defect_entry.anchor_atom_coords,
                                         np.array([0.5, 0.5, 0.5]))


def test_msonable(defect_entry):
    assert_msonable(defect_entry)


def test_json_roundtrip(defect_entry, tmpdir):
    assert_json_roundtrip(defect_entry, tmpdir)


def test_defect_entry(defect_entry, tmpdir):
    tmpdir.chdir()
    defect_entry.to_prior_info(filename="a.yaml")
    expected = """charge: 1
"""
    assert Path("a.yaml").read_text() == expected

"""
TODO
- Store basic info for a single input defect initial structure.
- create prior_info
- from structure
DONE
"""