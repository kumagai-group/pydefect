# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from monty.serialization import loadfn
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator, SiteDiff
from pydefect.tests.helpers.assertion import assert_msonable
from pymatgen import Structure, IStructure, Lattice


@pytest.fixture
def structure_comparator():
    cu2o_perfect = IStructure(Lattice.cubic(5),
                              species=["Cu"] * 4 + ["O"] * 2,
                              coords=[[0.25, 0.25, 0.25],  # removed
                                      [0.25, 0.74, 0.74],  # removed
                                      [0.75, 0.75, 0.25],
                                      [0.75, 0.25, 0.75],
                                      [0, 0, 0],
                                      [0.5, 0.5, 0.5]])
    cu2o_defect = IStructure(Lattice.cubic(5),
                             species=["Cu"] * 3 + ["O"] * 2 + ["H"],
                             coords=[[0.25, 0.5, 0.5],   # inserted
                                     [0.76, 0.73, 0.24],
                                     [0.75, 0.25, 0.73],
                                     [0.05, 0.95, 0],
                                     [0.5, 0.5, 0.5],
                                     [0.25]*3])   # inserted

    return DefectStructureComparator(defect_structure=cu2o_defect,
                                     perfect_structure=cu2o_perfect)


@pytest.fixture
def structure_comparator_periodic_issue():
    cu2o_perfect = IStructure(Lattice.cubic(5),
                              species=["Cu"] * 4 + ["O"] * 2,
                              coords=[[0.25, 0.25, 0.25],
                                      [0.25, 0.75, 0.75],
                                      [0.75, 0.75, 0.25],
                                      [0.75, 0.25, 0.75],
                                      [0, 0, 0],
                                      [0.5, 0.5, 0.5]])
    # defect center is ([1.0, 1.0, 1.0] + [0.99, 0.99, 0.99]) / 2 = [0.995]*3
    cu2o_defect = IStructure(Lattice.cubic(5),
                             species=["Cu"] * 4 + ["O"] + ["H"],
                             coords=[[0.25, 0.25, 0.25],
                                     [0.25, 0.75, 0.75],
                                     [0.75, 0.75, 0.25],
                                     [0.75, 0.25, 0.75],
                                     [0.5, 0.5, 0.5],
                                     [0.99, 0.99, 0.99]])

    return DefectStructureComparator(defect_structure=cu2o_defect,
                                     perfect_structure=cu2o_perfect)


def test_atom_mapping_to_perfect(structure_comparator):
    assert structure_comparator.atom_mapping == {1: 2, 2: 3, 3: 4, 4: 5}
    assert structure_comparator.removed_indices == [0, 1]
    assert structure_comparator.inserted_indices == [0, 5]


def test_defect_structure_analyzer_defect_center(structure_comparator):
    actual = structure_comparator.defect_center_coord
    assert (actual == np.array([0.25, 0.435, 0.435])).all()


def test_defect_center_periodicity(structure_comparator_periodic_issue):
    actual = structure_comparator_periodic_issue.defect_center_coord
    assert (actual == np.array([0.995]*3)).all()


def test_neighboring_atom_indices(cubic_supercell):
    structure: Structure = cubic_supercell.copy()
    structure.pop(32)  # [0.25, 0, 0]
    structure.pop(0)  # [0, 0, 0]
    structure.append(species="Li", coords=[0.124, 0, 0])
    structure.to(filename="POSCAR")
    structure_analyzer = DefectStructureComparator(structure, cubic_supercell)
    actual = structure_analyzer.neighboring_atom_indices()
    expected = sorted([25, 15, 16, 23, 54, 47, 46, 56, 62])
    assert actual == expected


def test_actual_files(vasp_files):
    d = vasp_files / "KInO2_Va_O_2"
    calc_results = loadfn(d / "Va_O1_2/calc_results.json")
    perfect_calc_results = loadfn(d / "perfect_calc_results.json")
    structure_comparator = DefectStructureComparator(
        calc_results.structure, perfect_calc_results.structure)
    vac_atom_indices = list(range(192))
    vac_atom_indices.pop(96)
    expected = dict(zip(range(191), vac_atom_indices))
    assert structure_comparator.inserted_indices == []
    assert structure_comparator.removed_indices == [96]
    assert structure_comparator.atom_mapping == expected


def test_site_diff():
    site_diff = SiteDiff(removed=[("H", 0, (0.0, 0.0, 0.0))],
                         inserted=[("H", 0, (0.0, 0.0, 0.0))],
                         mapping=[None, 0, 1, 2])
    assert_msonable(site_diff)


# def test_make_defect_type():
#     lattice = Lattice.cubic(10)
#     assert make_defect_type([("")], [], lattice) == DefectType.vacancy
#     assert make_defect_type([], [1], lattice) == DefectType.interstitial
#     assert make_defect_type([3], [1], lattice) == DefectType.interstitial

