# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.input_maker.defect_set import DefectSet, DefectSetMaker, charge_set
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


def test_charge_set():
    assert charge_set(3) == [-1, 0, 1, 2, 3]
    assert charge_set(2) == [0, 1, 2]
    assert charge_set(-3) == [-3, -2, -1, 0, 1]
    assert charge_set(-2) == [-2, -1, 0]


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

#
# def test_defect_set_maker_dopants(supercell_info):
#     defect_maker = DefectSetMaker(supercell_info, dopants=["Li"])
#     assert defect_maker.defect_set.names == [
#         DefectName.from_str("Va_H1_-1"),
#         DefectName.from_str("Va_H1_0"),
#         DefectName.from_str("Va_H1_1"),
#         DefectName.from_str("Va_He1_0"),
#         DefectName.from_str("Li_H1_0"),
#         DefectName.from_str("Li_He1_-1"),
#         DefectName.from_str("Li_He1_0"),
#         DefectName.from_str("Li_He1_1"),
#     ]
#

"""
TODO
- Dopant

DONE
- Set Defect names
- Write yaml files
"""

