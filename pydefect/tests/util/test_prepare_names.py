# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.util.prepare_names import defect_mpl_name, defect_plotly_name, \
    defect_plotly_full_name


def test_defect_mpl_name():
    assert defect_mpl_name(name="Va_O1") == "$V_{{\\rm O}1}$"
    assert defect_mpl_name(name="Mg_i1") == "${\\rm Mg}_{i1}$"


def test_defect_plotly_name():
    assert defect_plotly_name(name="Va_O1") == "<i>V</i><sub>O1</sub>"


def test_defect_plotly_fullname():
    actual = defect_plotly_full_name(fullname="Va_O1_1")
    assert actual == "<i>V</i><sub>O1</sub><sup>1</sup>"


