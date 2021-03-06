# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy

import numpy as np
import pytest
from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import \
    Displacement, DefectStructureInfo, fold_coords_in_structure, calc_drift, \
    make_defect_structure_info, symmetry_relation, SymmRelation, \
    unique_point_group, elem_indices_coords, DefectType, judge_defect_type, \
    calc_displacements
from pymatgen.core import Structure, Lattice
from vise.tests.helpers.assertion import assert_dataclass_almost_equal, \
    assert_msonable, assert_json_roundtrip
from vise.util.enum import ExtendedEnum
from vise.util.structure_symmetrizer import num_sym_op


@pytest.fixture
def structures():
    perfect = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Li", "U"],
        coords=[[0.25, 0.25, 0.25],
                [0.75, 0.75, 0.75],
                [0.5, 0.5, 0.5],
                [0, 0, 0]])

    initial = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.25, 0.25, 0.25],
                [0.75, 0.75, 0.75],
                [0.5, 0.5, 0.5],
                [0.0, 0.0, 0.0]])

    final = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.27, 0.25, 0.2501],
                [0.76, 0.75, 0.7501],
                [0.5, 0.5, 0.501],
                [0.0, 0.0, 0.0001]])
    return perfect, initial, final


@pytest.fixture
def displacements():
    return [None,
            Displacement(specie="He",
                         original_pos=(0.75, 0.75, 0.75),
                         final_pos=(0.76, 0.75, 0.75),
                         distance_from_defect=6.465437190360757,
                         disp_vector=(0.1, 0.0, 0.0),
                         displace_distance=0.1,
                         angle=124.9),
            Displacement(specie="Be",
                         original_pos=(0.5, 0.5, 0.5),
                         final_pos=(0.5, 0.5, 0.5009),
                         distance_from_defect=2.1355627039494767,
                         disp_vector=(0.0, 0.0, 0.009),
                         displace_distance=0.009,
                         angle=125.8),
            Displacement(specie="U",
                         original_pos=(0, 0, 0),
                         final_pos=(0.0, 0.0, 0.0),
                         distance_from_defect=6.525191802736529,
                         disp_vector=(0.0, 0.0, 0.0),
                         displace_distance=0.0,
                         angle=None)]


@pytest.fixture
def def_str_info(displacements):
    site_diff = SiteDiff(removed=[(0, 'H', (0.25, 0.25, 0.25))],
                         inserted=[(0, 'H', (0.27, 0.25, 0.25))],
                         removed_by_sub=[(2, 'Li', (0.5, 0.5, 0.5))],
                         inserted_by_sub=[(2, 'Be', (0.5, 0.5, 0.5009))])
    site_diff_from_init = SiteDiff(removed=[(0, 'H', (0.25, 0.25, 0.25))],
                                   inserted=[(0, 'H', (0.27, 0.25, 0.25))],
                                   removed_by_sub=[(2, 'Li', (0.5, 0.5, 0.5))],
                                   inserted_by_sub=[(2, 'Be', (0.5, 0.5, 0.5))])
    return DefectStructureInfo(initial_site_sym="3m",
                               final_site_sym="m",
                               site_diff=site_diff,
                               site_diff_from_initial=site_diff_from_init,
                               symprec=0.1,
                               dist_tol=0.2,
                               anchor_atom_idx=3,
                               drift_vector=(0.0, 0.0, 0.0001),
                               drift_dist=0.001,
                               center=(0.38, 0.375, 0.375175),
                               displacements=displacements)


def test_fold_coords():
    actual = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                       species=["H"],
                       coords=[[0.5, 0, 0.49]])
    fold_coords_in_structure(actual, [0, 0, 0])
    expected = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                         species=["H"],
                         coords=[[-0.5, 0, 0.49]])
    assert actual == expected


def test_calc_drift_dist(structures):
    perfect, _, final = structures
    anchor_atom_idx, drift_dist, vector_in_frac = \
        calc_drift(perfect, final, center=[0.4, 0.4, 0.4], d_to_p=[None, 1, 2, 3])
    assert anchor_atom_idx == 3
    assert drift_dist == 0.001
    np.testing.assert_almost_equal(vector_in_frac, np.array([0.0, 0.0, 0.0001]))


def test_calc_displacements():
    perf = Structure(Lattice.cubic(10), species=["H"], coords=[[0, 0, -0.01]])
    defect = Structure(Lattice.cubic(10), species=["H"], coords=[[0, 0, 10.98]])
    actual = calc_displacements(perf, defect, [0, 0, 0.98], d_to_p=[0])
    expected = [Displacement(specie="H",
                             original_pos=(0.0, 0.0, 0.99),
                             final_pos=(0.0, 0.0, 0.98),
                             distance_from_defect=0.1,
                             disp_vector=(0.0, 0.0, -0.1),
                             displace_distance=0.1,
                             angle=0.0)]
    assert_dataclass_almost_equal(actual[0], expected[0])


def test_elem_indices_coords(structures):
    actual = elem_indices_coords(structures[0], indices=[1, 2])
    expected = [("He", 1, (0.75, 0.75, 0.75)), ("Li", 2, (0.5, 0.5, 0.5))]
    assert actual == expected


def test_defect_type():
    assert_msonable(DefectType.vacancy)
    assert issubclass(DefectType, ExtendedEnum)


def test_judge_defect_type():
    """The structures remain as they could be used in the future"""
    perfect = Structure(Lattice.cubic(10), species=["H"]*4,
                        coords=[[0.0, 0.0, 0.0], [0.0, 0.5, 0.5],
                                [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    defect = Structure(Lattice.cubic(10), species=["H"]*3,
                       coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.5],
                               [0.5, 0.5, 0.0]])
    site_diff = SiteDiff(removed=[(1, "H", (0.0, 0.5, 0.5)),
                                  (2, "H", (0.5, 0.0, 0.5))],
                         inserted=[(1, "H", (0.25, 0.25, 0.5))],
                         removed_by_sub=[], inserted_by_sub=[])
    assert judge_defect_type(site_diff) == DefectType.vacancy_split

    defect_2 = Structure(Lattice.cubic(10), species=["H"]*5,
                         coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.5],
                                 [0.75, 0.25, 0.5], [0.25, 0.75, 0.5],
                                 [0.5, 0.5, 0.0]])
    site_diff_2 = SiteDiff(removed=[(1, "H", (0.0, 0.5, 0.5)),
                                    (2, "H", (0.5, 0.0, 0.5))],
                           inserted=[(1, "H", (0.25, 0.25, 0.5)),
                                     (2, "H", (0.75, 0.25, 0.5)),
                                     (3, "H", (0.25, 0.75, 0.5))],
                           removed_by_sub=[], inserted_by_sub=[])
    assert judge_defect_type(site_diff_2) == DefectType.interstitial_split

    defect_3 = Structure(Lattice.cubic(10), species=["H", "He", "H"],
                         coords=[[0.0, 0.0, 0.0], [0.25, 0.25, 0.5],
                                 [0.5, 0.5, 0.0]])
    site_diff_3 = SiteDiff(removed=[(1, "H", (0.0, 0.5, 0.5)),
                                    (2, "H", (0.5, 0.0, 0.5))],
                           inserted=[(1, "He", (0.25, 0.25, 0.5))],
                           removed_by_sub=[], inserted_by_sub=[])

    assert judge_defect_type(site_diff_3) == DefectType.unknown


def test_str_info_msonable(def_str_info):
    assert_msonable(def_str_info)


def test_calc_results_json_roundtrip(def_str_info, tmpdir):
    assert_json_roundtrip(def_str_info, tmpdir)


def test_symmetry_relation():
    for point_group in num_sym_op:
        pg = unique_point_group(point_group)
        assert symmetry_relation(pg, pg) == SymmRelation.same
    assert symmetry_relation("3m", "m") == SymmRelation.subgroup
    assert symmetry_relation("4", "2") == SymmRelation.subgroup
    assert symmetry_relation("-43m", "mm2") == SymmRelation.subgroup
    assert symmetry_relation("4/m", "2") == SymmRelation.subgroup


def test_make_defect_structure_info(structures, def_str_info):
    perfect, initial, final = structures
    actual = make_defect_structure_info(
         perfect, initial, final, dist_tol=0.2, symprec=0.1)
    assert_dataclass_almost_equal(actual, def_str_info)


def test_make_defect_structure_info_w_symms(structures, def_str_info):
    perfect, initial, final = structures
    actual = make_defect_structure_info(
        perfect, initial, final, dist_tol=0.2,
        init_site_sym="3m", final_site_sym="m")
    expected = deepcopy(def_str_info)
    expected.symprec = None
    assert_dataclass_almost_equal(actual, expected)


def test_make_def_str_info_symm_rel(structures, def_str_info):
    assert def_str_info.symm_relation == SymmRelation.subgroup


def test_make_def_str_info_same_config_from_init(structures, def_str_info):
    assert def_str_info.same_config_from_init is False


def test_def_str_info_defect_type(def_str_info):
    assert def_str_info.defect_type == DefectType.unknown


def test_repr(def_str_info):
    print(def_str_info)


