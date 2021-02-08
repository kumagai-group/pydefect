# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pydefect.analyzer.defect_structure_info import \
    Displacement, DefectStructureInfo, fold_coords, calc_drift, \
    calc_displacements, make_defect_structure_info, symmetrize_defect_structure, \
    symmetry_relation, SymmRelation, unique_poing_group, elem_and_indices
from pydefect.defaults import defaults
from pydefect.tests.helpers.assertion import assert_msonable
from pymatgen import Structure, Lattice
from vise.util.structure_symmetrizer import StructureSymmetrizer
from vise.util.structure_symmetrizer import num_sym_op


@pytest.fixture
def structures():
    perfect = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Li", "U"],
        coords=[[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0, 0, 0]])

    initial = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.55, 0.5, 0.5], [0.0, 0.0, 0.0]])

    final = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.35, 0.35, 0.35], [0.74, 0.74, 0.74], [0.5, 0.5, 0.5], [0.0, 0.0, 0.0001]])
    return perfect, initial, final


@pytest.fixture
def displacements():
    x = np.sqrt((0.75-0.4)**2*3) * 10
    return [None,
            Displacement(specie="He",
                         original_pos=(0.75, 0.75, 0.75),
                         final_pos=(0.74, 0.74, 0.74),
                         distance_from_defect=round(x, 14),
                         displace_distance=0.17320508075688712,
                         angle=0.0),
            None,
            Displacement(specie="U",
                         original_pos=(0, 0, 0),
                         final_pos=(0.0, 0.0, 0.0001),
                         distance_from_defect=np.sqrt(3) * 10 * 0.4,
                         displace_distance=0.001,
                         angle=54.7)]


@pytest.fixture
def def_str_info(structures, displacements):
    perfect, initial, final = structures
    return DefectStructureInfo(initial_point_group="m-3m",
                               final_point_group="m-3m",
                               initial_structure=initial,
                               final_structure=final,
                               perfect_structure=perfect,
                               fd_to_p=[None, 1, None, 3],
#                               fd_to_id=[None, 1, 2, 3],
                               drift_dist=0.001,
#                               initial_defect_type=DefectType.substituted,
                               initial_vacancies=[("Li", 2)],
                               initial_interstitials=[("Be", 2)],
#                               final_defect_type=DefectType.complex,
                               final_vacancies=[("H", 0), ("Li", 2)],
                               final_interstitials=[("H", 0), ("Be", 2)],
                               defect_center_coord=(0.4, 0.4, 0.4),
                               displacements=displacements,
                               neighboring_atom_indices=[0, 1, 2])


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
    drift_dist, vector_in_frac = calc_drift(perfect, final,
                                            center=[0.4, 0.4, 0.4],
                                            d_to_p=[None, 1, 2, 3])
    assert drift_dist == 0.001
    np.testing.assert_almost_equal(vector_in_frac, np.array([0.0, 0.0, 0.0001]))


def test_calc_displacements(structures, displacements):
    perfect, _, final = structures
    actual = calc_displacements(perfect, final,
                                center=[0.4, 0.4, 0.4],
                                d_to_p=[None, 1, None, 3])
    expected = displacements
    assert actual == expected


def test_elem_and_indices(structures):
    actual = elem_and_indices(structures[0], indices=[1, 2])
    expected = [("He", 1), ("Li", 2)]
    assert actual == expected


def test_str_info_msonable(def_str_info):
    assert_msonable(def_str_info)


def test_symmetry_relation():
    for point_group in num_sym_op:
        pg = unique_poing_group(point_group)
        assert symmetry_relation(pg, pg) == SymmRelation.same
    assert symmetry_relation("3m", "m") == SymmRelation.subgroup
    assert symmetry_relation("4", "2") == SymmRelation.subgroup
    assert symmetry_relation("-43m", "mm2") == SymmRelation.subgroup
    assert symmetry_relation("4/m", "2") == SymmRelation.subgroup


def test_make_defect_structure_info(structures, def_str_info, mocker):
    perfect, initial, final = structures
    p_calc_results = mocker.stub(name='perfect_calc_results')
    defect_entry = mocker.stub(name='defect_entry')
    d_calc_results = mocker.stub(name='calc_results')

    p_calc_results.structure = perfect
    defect_entry.structure = initial
    defect_entry.site_symmetry = "m-3m"
    d_calc_results.structure = final
    d_calc_results.site_symmetry = "m-3m"
    actual = make_defect_structure_info(
        p_calc_results, defect_entry, d_calc_results, neighbor_cutoff_factor=3)
    assert actual.__repr__() == def_str_info.__repr__()


def test_repr(def_str_info):
    print(def_str_info)


