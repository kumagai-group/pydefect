# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from pydefect.input_maker.defect_entry import DefectEntry, make_defect_entry, \
    PerturbedSite
from pymatgen.core import Lattice, IStructure

# "H" at [0.0, 0.0, 0.0] is removed here.
from vise.tests.helpers.assertion import assert_msonable, assert_json_roundtrip

perf_coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
               [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.5]]

coords = deepcopy(perf_coords)
coords.pop(0)

perturbed_coords = \
    [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
     [0.0, 0.0, 0.501], [0.0, 0.501, 0.0], [0.501, 0.0, 0.0], [0.5, 0.5, 0.5]]

perfect = IStructure(
    Lattice.cubic(10.0), ["H"] * 4 + ["He"] * 4, perf_coords)
defect_structure = IStructure(
    Lattice.cubic(10.0), ["H"] * 3 + ["He"] * 4, coords)
perturbed_defect = IStructure(
    Lattice.cubic(10.0), ["H"] * 3 + ["He"] * 4, perturbed_coords)


@pytest.fixture
def defect_entry():
    perturbed_sites = (
        PerturbedSite("He", 5.0, (0.0, 0.0, 0.5), (0.0, 0.0, 0.501), 0.1),
        PerturbedSite("He", 5.0, (0.0, 0.5, 0.0), (0.0, 0.501, 0.0), 0.1),
        PerturbedSite("He", 5.0, (0.5, 0.0, 0.0), (0.501, 0.0, 0.0), 0.1))
    return DefectEntry(name="Va_O1",
                       charge=1,
                       structure=defect_structure,
                       site_symmetry="m-3m",
                       defect_center=(0.0, 0.0, 0.0),
                       perturbed_structure=perturbed_defect,
                       perturbed_sites=perturbed_sites,
                       perturbed_site_symmetry="1")


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


def test_defect_entry_str(defect_entry, tmpdir):
    actual = defect_entry.__str__()
    expected = """ -- defect entry info
name: Va_O1_1
site symmetry: m-3m
defect center: ( 0.000,  0.000,  0.000)
perturbed sites:
elem dist   initial_coords             perturbed_coords         displacement
  He 5.00 ( 0.000,  0.000,  0.500) -> ( 0.000,  0.000,  0.501)    0.10
  He 5.00 ( 0.000,  0.500,  0.000) -> ( 0.000,  0.501,  0.000)    0.10
  He 5.00 ( 0.500,  0.000,  0.000) -> ( 0.501,  0.000,  0.000)    0.10"""
    assert actual == expected


def test_make_defect_entry(defect_entry):
    relaxed_coords = \
        [[0.25, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
         [0.0, 0.0, 0.51], [0.0, 0.51, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.5]]
    relaxed_defect = IStructure(
        Lattice.cubic(10.0), ["Li"] + ["H"] * 3 + ["He"] * 4, relaxed_coords)

    unrelaxed_coords = \
        [[0.25, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
         [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.5]]

    unrelaxed_defect = IStructure(
        Lattice.cubic(10.0), ["Li"] + ["H"] * 3 + ["He"] * 4, unrelaxed_coords)

    actual = make_defect_entry(name="Va_O1",
                               charge=1,
                               perfect_structure=perfect,
                               defect_structure=relaxed_defect)

    expected = DefectEntry(name="Va_O1",
                           charge=1,
                           structure=unrelaxed_defect,
                           perturbed_structure=None,
                           site_symmetry="4mm",
                           defect_center=(0.125, 0.0, 0.0))

    assert actual == expected
