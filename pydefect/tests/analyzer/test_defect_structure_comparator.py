# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from monty.serialization import loadfn
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator, SiteDiff
from pymatgen.core import Structure, IStructure, Lattice
from vise.tests.helpers.assertion import assert_msonable


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


@pytest.fixture
def site_diff():
    return SiteDiff(removed=[(1, "Cu", (0.25, 0.74, 0.74))],
                    inserted=[(0,  "Cu", (0.25, 0.5, 0.5))],
                    removed_by_sub=[(0,  "Cu", (0.25, 0.25, 0.25))],
                    inserted_by_sub=[(5, "H", (0.25, 0.25, 0.25))])


def test_site_diff_msonable(site_diff):
    assert_msonable(site_diff)


def test_site_diff_is_no_diff():
    site_diff = SiteDiff(removed=[], inserted=[],
                         removed_by_sub=[], inserted_by_sub=[])
    assert site_diff.is_no_diff is True
    site_diff = SiteDiff(removed=[(1, "H", (0.0, 0.0, 0.0))], inserted=[],
                         removed_by_sub=[], inserted_by_sub=[])
    assert site_diff.is_no_diff is False


def test_site_diff():
    site_diff_vac = SiteDiff(removed=[(1, "H", (0.0, 0.0, 0.0))], inserted=[],
                             removed_by_sub=[], inserted_by_sub=[])
    site_diff_int = SiteDiff(removed=[], inserted=[(1, "H", (0.0, 0.0, 0.0))],
                             removed_by_sub=[], inserted_by_sub=[])
    site_diff_sub = SiteDiff(removed=[], inserted=[],
                             removed_by_sub=[(1, "H", (0.0, 0.0, 0.0))],
                             inserted_by_sub=[(1, "He", (0.0, 0.0, 0.0))])
    assert site_diff_vac.is_complex_defect is False
    assert site_diff_sub.is_complex_defect is False

    assert site_diff_vac.is_vacancy
    assert site_diff_int.is_vacancy is False
    assert site_diff_sub.is_vacancy is False

    assert site_diff_vac.is_interstitial is False
    assert site_diff_int.is_interstitial
    assert site_diff_sub.is_interstitial is False

    assert site_diff_vac.is_substituted is False
    assert site_diff_int.is_substituted is False
    assert site_diff_sub.is_substituted


def test_make_site_diff(structure_comparator, site_diff):
    assert structure_comparator.make_site_diff() == site_diff


def test_make_site_diff_wo_diff():
    s1 = IStructure(Lattice.cubic(10), species=["H"], coords=[[0]*3])
    s2 = s1.copy()
    structure_comparator = DefectStructureComparator(s1, s2)
    site_diff = structure_comparator.make_site_diff()
    assert site_diff == SiteDiff([], [], [], [])
    assert site_diff.is_no_diff

    s1 = IStructure(Lattice.cubic(10), species=["H", "He"], coords=[[0]*3, [0.5]*3])
    s2 = IStructure(Lattice.cubic(10), species=["H"], coords=[[0]*3])
    structure_comparator = DefectStructureComparator(defect_structure=s2,
                                                     perfect_structure=s1)
    site_diff = structure_comparator.make_site_diff()
    assert site_diff == SiteDiff([(1, "He", (0.5, 0.5, 0.5))], [], [], [])
    assert site_diff.is_no_diff is False



# def test_make_defect_type():
#     lattice = Lattice.cubic(10)
#     assert make_defect_type([("")], [], lattice) == DefectType.vacancy
#     assert make_defect_type([], [1], lattice) == DefectType.interstitial
#     assert make_defect_type([3], [1], lattice) == DefectType.interstitial

