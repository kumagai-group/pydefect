# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_set_maker import charge_set
from pydefect.input_maker.defect import SimpleDefect


simple_defects = {SimpleDefect(None, "O1", [1]), SimpleDefect("N", "O1", [0, 1])}


@pytest.fixture
def defect_set():
    return DefectSet(defects=simple_defects)


def test_defect_set(defect_set):
    assert SimpleDefect(None, "O1", [1]) in defect_set
    assert len(defect_set) == 2
    for defect in defect_set:
        defect


def test_to_yaml(defect_set, tmpdir):
    tmpdir.chdir()
    defect_set.to_yaml("tmp.yaml")
    expected = """N_O1: [0, 1]
Va_O1: [1]
"""
    assert tmpdir.join("tmp.yaml").read() == expected


def test_from_file(defect_set, tmpdir):
    tmpdir.chdir()
    tmpdir.join("tmp.yaml").write("""N_O1: [0, 1]
Va_O1: [1] 
""")
    defect_set = DefectSet.from_yaml("tmp.yaml")
    assert isinstance(defect_set, DefectSet)
    assert defect_set == simple_defects


def test_charge_set():
    assert charge_set(3) == [-1, 0, 1, 2, 3]
    assert charge_set(2) == [0, 1, 2]
    assert charge_set(-3) == [-3, -2, -1, 0, 1]
    assert charge_set(-2) == [-2, -1, 0]


"""
TODO
- Dopant

DONE
- Set Defect names
- Write yaml files
"""
