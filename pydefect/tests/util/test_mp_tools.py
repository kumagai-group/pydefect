# -*- coding: utf-8 -*-

from pydefect.util.mp_tools import MpQuery, chemsys


def test_mp_actual_query():
    element_list = ["Mg", "O"]
    mp_query = MpQuery(element_list=element_list,
                       properties=["formula_pretty", "structure"])
    print([m.formula_pretty for m in mp_query.materials])


def test_chemsys():
    assert (chemsys(["Mg", "Al", "O"]) ==
            ["Mg", "Al", "O", 'Mg-Al', 'Mg-O', 'Al-O', 'Mg-Al-O'])
