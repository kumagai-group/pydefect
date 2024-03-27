# -*- coding: utf-8 -*-

from pydefect.util.mp_tools import MpQuery


def test_mp_actual_query():
    element_list = ["Mg", "O"]
    mp_query = MpQuery(element_list=element_list,
                       properties=["formula_pretty", "structure"])
    print(mp_query.materials)


