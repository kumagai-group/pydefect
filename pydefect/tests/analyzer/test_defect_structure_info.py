# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy

from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import \
    symmetry_relation, SymmRelation, \
    DefectType, judge_defect_type, unique_point_group
from pymatgen.core import Structure, Lattice
from vise.tests.helpers.assertion import assert_msonable, assert_json_roundtrip
from vise.util.enum import ExtendedEnum
from vise.util.structure_symmetrizer import num_sym_op


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


def test_defect_structure_info_json_roundtrip(def_str_info, tmpdir):
    assert_json_roundtrip(def_str_info, tmpdir)


def test_symmetry_relation():
    for point_group in num_sym_op:
        pg = unique_point_group(point_group)
        assert symmetry_relation(pg, pg) == SymmRelation.same
    assert symmetry_relation("3m", "m") == SymmRelation.subgroup
    assert symmetry_relation("4", "2") == SymmRelation.subgroup
    assert symmetry_relation("-43m", "mm2") == SymmRelation.subgroup
    assert symmetry_relation("4/m", "2") == SymmRelation.subgroup


def test_make_def_str_info_symm_rel(def_str_info):
    assert def_str_info.symm_relation == SymmRelation.subgroup


def test_make_def_str_info_same_config_from_init(def_str_info):
    assert def_str_info.same_config_from_init is False


def test_def_str_info_defect_type(def_str_info):
    assert def_str_info.defect_type == DefectType.unknown


def test_repr(def_str_info):
    x = deepcopy(def_str_info)
    x.site_diff_from_initial.inserted = []
    actual = x.__str__()
    expected = """ -- defect structure info
Defect type: unknown
Site symmetry: 3m -> m (subgroup)
Has same configuration from initial structure: False
Drift distance: 0.001
Defect center: ( 0.380,  0.375,  0.375)
Removed atoms:
0  H  ( 0.250,  0.250,  0.250)

Added atoms:
0  H  ( 0.270,  0.250,  0.250)

Removed atoms from initial structure:
0  H  ( 0.250,  0.250,  0.250)

Displacements
Elem  Dist  Displace  Angle  Index  Initial site                  Final site
Be    2.14  0.01      130    2      ( 0.500,  0.500,  0.500)  ->  ( 0.500,  0.500,  0.501)"""
    assert actual == expected


