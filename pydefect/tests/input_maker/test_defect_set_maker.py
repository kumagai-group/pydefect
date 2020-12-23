# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_set_maker import DefectSetMaker


def test_defect_set_maker(supercell_info):
    defect_maker = DefectSetMaker(supercell_info)
    expected = {SimpleDefect(None, "H1", [-1, 0, 1]),
                SimpleDefect(None, "He1", [0]),
                SimpleDefect("H", "i1", [-1, 0, 1]),
                SimpleDefect("He", "i1", [0])}
    assert defect_maker.defect_set == expected


def test_defect_set_maker_oxi_state(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, {"H": 2})
    assert defect_maker.defect_set == {SimpleDefect(None, "H1", [-2, -1, 0]),
                                                           SimpleDefect(None, "He1", [0]),
                                                           SimpleDefect("H", "i1", [0, 1, 2]),
                                                           SimpleDefect("He", "i1", [0])}


def test_defect_set_maker_dopants(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, dopants=["Li"])
    assert set(defect_maker.defect_set) == {
        SimpleDefect(None, "H1", [-1, 0, 1]),
        SimpleDefect(None, "He1", [0]),
        SimpleDefect("Li", "H1", [0]),
        SimpleDefect("H", "i1", [-1, 0, 1]),
        SimpleDefect("He", "i1", [0]),
        SimpleDefect("Li", "i1", [-1, 0, 1]),
    }


def test_defect_set_maker_electroneg_diff(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, dopants=["Li"],
                                  ele_neg_diff=2.2 - 0.98 - 1e-5)
    assert set(defect_maker.defect_set) == {
        SimpleDefect(None, "H1", [-1, 0, 1]),
        SimpleDefect(None, "He1", [0]),
        SimpleDefect("H", "i1", [-1, 0, 1]),
        SimpleDefect("He", "i1", [0]),
        SimpleDefect("Li", "i1", [-1, 0, 1]),
    }


def test_defect_set_maker_screening(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, keywords=["Va_H1"])
    expected = {SimpleDefect(None, "H1", [-1, 0, 1])}
    assert defect_maker.defect_set == expected


