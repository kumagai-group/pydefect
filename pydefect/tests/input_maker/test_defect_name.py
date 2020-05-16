# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest


from pydefect.input_maker.defect_name import DefectName
from pydefect.tests.helpers.assertion import assert_msonable


def test_defect_name_msonable():
    defect_name = DefectName("Va", "O1", 1)
    assert_msonable(defect_name)


def test_defect_name():
    defect_name = DefectName("Va", "O1", 1)
    assert defect_name.in_atom == "Va"
    assert defect_name.out_atom == "O1"
    assert defect_name.charge == 1
    assert str(defect_name) == "Va_O1_1"


def test_defect_name_hashable():
    defect_name = DefectName("Va", "O1", 1)
    d = {defect_name: 1}


def test_defect_name_eq():
    defect_name_1 = DefectName("Va", "O1", 1)
    defect_name_2 = DefectName("Va", "O1", 1)
    assert defect_name_1 == defect_name_2


def test_defect_name_from_str():
    actual = DefectName.from_str("Va_O1_1")
    expected = DefectName("Va", "O1", 1)
    assert actual == expected


def test_defect_name_specie():
    name = DefectName.from_str("Va_O1_1")
    assert name.name_wo_charge == "Va_O1"


"""
TODO
- 

DONE
"""