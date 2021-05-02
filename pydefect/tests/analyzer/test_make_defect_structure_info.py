# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import numpy as np
import pytest
from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import Displacement, \
    DefectStructureInfo
from pydefect.analyzer.make_defect_structure_info import \
    MakeDefectStructureInfo
from pydefect.cli.make_defect_vesta_file import fold_coords_in_structure

from pymatgen import Structure, Lattice
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


def test_fold_coords():
    actual = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                       species=["H"],
                       coords=[[0.5, 0, 0.49]])
    fold_coords_in_structure(actual, [0, 0, 0])
    expected = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                         species=["H"],
                         coords=[[-0.5, 0, 0.49]])
    assert actual == expected


def test_make_defect_structure_info(structures, def_str_info):
    perf, initial, final = structures
    info = MakeDefectStructureInfo(perf, initial, final,  dist_tol=0.2,
                                   symprec=0.1,
                                   neighbor_cutoff_factor=1.2)
    actual = info.defect_structure_info
    print(actual)
    assert_dataclass_almost_equal(actual, def_str_info)


def test_make_defect_structure_info2():
    perf = Structure(Lattice.cubic(1), ["H"], [[0., 0., 0.]])
    init = Structure(Lattice.cubic(1), ["H"]*2, [[0., 0., 0.], [0.5, 0.0, 0.0]])
    fin = Structure(Lattice.cubic(1), ["H"]*2, [[0., 0., 0.], [0.4, 0.0, 0.0]])

    info = MakeDefectStructureInfo(perf, init, fin,  dist_tol=0.05,
                                   symprec=0.1,
                                   neighbor_cutoff_factor=1.2)
    actual = info.defect_structure_info
    print(actual)

