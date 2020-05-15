# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_name import DefectName


names = [DefectName.from_str("Va_O1_1"), DefectName.from_str("N_O1_0")]


@pytest.fixture
def defect_set():
    return DefectSet(names=names)


def test_set_names(defect_set):
    assert defect_set.names[0] == names[0]


def test_to_file(defect_set, tmpdir):
    tmpdir.chdir()
    defect_set.to_file("tmp")
    expected = """Va_O1_1
N_O1_0"""
    assert tmpdir.join("tmp").read() == expected


def test_from_file(defect_set, tmpdir):
    tmpdir.chdir()
    tmpdir.join("tmp").write("""Va_O1_1
N_O1_0""")
    defect_set = DefectSet.from_file("tmp")
    assert isinstance(defect_set, DefectSet)
    assert defect_set.names == names

"""
TODO
- Write yaml files

DONE
- Set Defect names
"""