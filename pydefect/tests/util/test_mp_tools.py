# -*- coding: utf-8 -*-

# from vise.util.mp_tools import make_poscars_from_mp
# .testing import ViseTest

from pydefect.util.mp_tools import elements, MpQuery


def test_elements():
    assert len(elements) == 118
    assert elements[0] == "H"
    assert elements[-1] == "Og"


def test_mp_actual_query():
    element_list = ["Mg", "O"]
    mp_query = MpQuery(element_list=element_list,
                       properties=["full_formula", "structure"])
    print(mp_query.materials)


def test_mp_query(mocker):
    mock = mocker.patch("pydefect.util.mp_tools.MPRester")
    expected = [{'full_formula': 'Mg9'}, {'full_formula': 'Mg1O1'}, {'full_formula': 'O8'}, {'full_formula': 'O2'}, {'full_formula': 'O2'}, {'full_formula': 'O8'}]
    mock.return_value.__enter__.return_value.query.return_value = expected
    element_list = ["Mg", "O"]
    excluded = list(set(elements) - set(element_list))
    mp_query = MpQuery(element_list=element_list,
                       properties=["full_formula"])
    mock.return_value.__enter__.return_value.query.assert_called_with(
        criteria={"elements": {"$in": element_list, "$nin": excluded},
                  "e_above_hull": {"$lte": 1e-5}},
        properties=["full_formula"])
    assert mp_query.materials == expected


"""
TODO
- Get materials' formula containing given elements via MPRester.

DONE
- create list of elements. 
"""