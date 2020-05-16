# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.input_maker.defect_name import DefectName
from pydefect.input_maker.defect_set_maker import DefectSetMaker


def test_defect_set_maker(supercell_info):
    defect_maker = DefectSetMaker(supercell_info)
    assert defect_maker.defect_set.names == [DefectName.from_str("Va_H1_-1"),
                                             DefectName.from_str("Va_H1_0"),
                                             DefectName.from_str("Va_H1_1"),
                                             DefectName.from_str("Va_He1_0")]


def test_defect_set_maker_oxi_state(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, {"H": 0})
    assert defect_maker.defect_set.names == [DefectName.from_str("Va_H1_0"),
                                             DefectName.from_str("Va_He1_0")]


def test_defect_set_maker_dopants(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, dopants=["Li"])
    assert set(defect_maker.defect_set.names) == {
        DefectName.from_str("Va_H1_-1"),
        DefectName.from_str("Va_H1_0"),
        DefectName.from_str("Va_H1_1"),
        DefectName.from_str("Va_He1_0"),
        DefectName.from_str("Li_H1_0"),
    }


def test_defect_set_maker_electroneg_diff(supercell_info):
    defect_maker = DefectSetMaker(supercell_info, dopants=["Li"],
                                  ele_neg_diff=2.2 - 0.98 - 1e-5)
    assert set(defect_maker.defect_set.names) == {
        DefectName.from_str("Va_H1_-1"),
        DefectName.from_str("Va_H1_0"),
        DefectName.from_str("Va_H1_1"),
        DefectName.from_str("Va_He1_0"),
    }