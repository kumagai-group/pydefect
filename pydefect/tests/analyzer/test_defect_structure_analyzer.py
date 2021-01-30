# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pydefect.analyzer.defect_structure_analyzer import \
    Displacement, DefectStructureInfo, fold_coords, calc_drift_dist, \
    calc_displacements, make_defect_structure_info
from pydefect.tests.helpers.assertion import assert_msonable
from pymatgen import Structure, Lattice


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
                         distance_from_defect=round(x, 3),
                         displace_distance=round(np.sqrt(3) * 10 * 0.01, 3),
                         angle=0.0,
                         annotation="inward"),
            None,
            Displacement(specie="U",
                         original_pos=(0, 0, 0),
                         final_pos=(0.0, 0.0, 0.0001),
                         distance_from_defect=round(np.sqrt(3) * 10 * 0.4, 3),
                         displace_distance=0.001,
                         angle=54.7,
                         annotation=None),
            ]


@pytest.fixture
def def_str_info(structures, displacements):
    perfect, initial, final = structures
    return DefectStructureInfo(initial_point_group="m-3m",
                               final_point_group="m-3m",
                               initial_structure=initial,
                               final_structure=final,
                               perfect_structure=perfect,
#                               fd_to_id=[None, 1, 2, 3],
                               drift_dist=0.001,
#                               initial_defect_type=DefectType.substituted,
                               initial_vacancies=[2],
                               initial_interstitials=[2],
#                               final_defect_type=DefectType.complex,
                               final_vacancies=[0, 2],
                               final_interstitials=[0, 2],
                               defect_center_coord=(0.4, 0.4, 0.4),
                               displacements=displacements)
#                               neighboring_atom_indices=[0, 1, 2])


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
    drift_dist = calc_drift_dist(perfect, final,
                                 center=[0.4, 0.4, 0.4],
                                 d_to_p=[None, 1, 2, 3])
    assert drift_dist == 0.001


def test_calc_displacements(structures, displacements):
    perfect, _, final = structures
    actual = calc_displacements(perfect, final,
                                center=[0.4, 0.4, 0.4],
                                d_to_p=[None, 1, None, 3])
    expected = displacements
    assert actual == expected


def test_str_info_msonable(def_str_info):
    assert_msonable(def_str_info)


def test_make_defect_structure_info(structures, def_str_info, mocker):
    perfect, initial, final = structures
    defect_entry = mocker.stub(name='defect_entry')
    calc_results = mocker.stub(name='calc_results')
    defect_entry.structure = initial
    defect_entry.site_symmetry = "m-3m"
    calc_results.structure = final
    calc_results.site_symmetry = "m-3m"
    actual = make_defect_structure_info(perfect, defect_entry, calc_results)
    assert actual == def_str_info

# def test_repr(structure_analyzer):
#     expected = """
# Site symmetry:
#   m-3m -> 3m
# Defect center:
#   [0.25, 0.435, 0.435]
# Transferred atoms:
#  - Cu [0.25, 0.25, 0.25]
#  - Cu [0.25, 0.74, 0.74]
#  + Cu [0.25, 0.5, 0.5]
#  +  H [0.25, 0.25, 0.25]
# Defect type:
#   Complex defect
# Sum of displacement distances:
#   3.45 A
# Displacements:
#   Cu 1 [0.25, 0.75, 0.75] -> 1 [0.25, 0.75, 0.75] ([0.25, 0.75, 0.75])  0.0 A  0°
#   Cu 2 [0.25, 0.75, 0.75] -> 1 [0.25, 0.75, 0.75] ([0.25, 0.75, 0.75])  0.0 A  0° inward
#     """
#     assert structure_analyzer.atom_mapping == {1: 2, 2: 3, 3: 4, 4: 5}

