# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import \
    Displacement, DefectStructureInfo, fold_coords, calc_drift, \
    make_defect_structure_info, symmetry_relation, SymmRelation, \
    unique_point_group, elem_indices_coords, defect_vesta_file
from pydefect.tests.helpers.assertion import assert_msonable
from pymatgen import Structure, Lattice
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
                         displace_distance=0.09999999999999964,
                         angle=124.9),
            Displacement(specie="Be",
                         original_pos=(0.5, 0.5, 0.5),
                         final_pos=(0.5, 0.5, 0.5009),
                         distance_from_defect=2.1355627039494767,
                         disp_vector=(0.0, 0.0, 0.009),
                         displace_distance=0.009000000000000341,
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
    site_diff = SiteDiff(removed={0: ('H', (0.25, 0.25, 0.25))},
                         inserted={0: ('H', (0.27, 0.25, 0.25))},
                         removed_by_sub={2: ('Li', (0.5, 0.5, 0.5))},
                         inserted_by_sub={2: ('Be', (0.5, 0.5, 0.5009))},
                         mapping={1: 1, 3: 3})
    site_diff_from_init = SiteDiff(removed={0: ('H', (0.25, 0.25, 0.25))},
                                   inserted={0: ('H', (0.27, 0.25, 0.25))},
                                   removed_by_sub={2: ('Li', (0.5, 0.5, 0.5))},
                                   inserted_by_sub={2: ('Be', (0.5, 0.5, 0.5))},
                                   mapping={1: 1, 3: 3})
    return DefectStructureInfo(initial_site_sym="3m",
                               final_site_sym="m",
                               site_diff=site_diff,
                               site_diff_from_initial=site_diff_from_init,
                               symprec=0.1,
                               dist_tol=0.2,
                               anchor_atom_idx=3,
                               drift_vector=(0.0, 0.0, 0.0001),
                               drift_dist=0.001,
                               center=(0.38, 0.375, 0.37517500000000004),
                               displacements=displacements)


def test_fold_coords():
    actual = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                       species=["H"],
                       coords=[[0.5, 0, 0.49]])
    fold_coords(actual, [0, 0, 0])
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

#
# def test_calc_displacements(structures, displacements):
#     perfect, _, final = structures
#     actual = calc_displacements(perfect, final,
#                                 center=[0.4, 0.4, 0.4],
#                                 d_to_p=[None, 1, None, 3])
#     expected = displacements
#     assert actual == expected


def test_elem_indices_coords(structures):
    actual = elem_indices_coords(structures[0], indices=[1, 2])
    expected = [("He", 1, (0.75, 0.75, 0.75)), ("Li", 2, (0.5, 0.5, 0.5))]
    assert actual == expected


def test_str_info_msonable(def_str_info):
    assert_msonable(def_str_info)


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
    print(actual.displacements)
#    assert actual.__repr__() == def_str_info.__repr__()


def test_make_def_str_info_symm_rel(structures, def_str_info):
    assert def_str_info.symm_relation == SymmRelation.subgroup


def test_make_def_str_info_same_config_from_init(structures, def_str_info):
    assert def_str_info.same_config_from_init is False


def test_repr(def_str_info):
    print(def_str_info)


