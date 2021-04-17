# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import pytest
from pydefect.util.prepare_names import defect_mpl_name, defect_plotly_name, \
    defect_plotly_full_name, typical_defect_name, prettify_names


def test_defect_mpl_name():
    assert defect_mpl_name(name="Va_O1") == "$V_{{\\rm O}1}$"
    assert defect_mpl_name(name="Mg_i1") == "${\\rm Mg}_{i1}$"


def test_defect_plotly_name():
    assert defect_plotly_name(name="Va_O1") == "<i>V</i><sub>O1</sub>"


def test_defect_plotly_fullname():
    actual = defect_plotly_full_name(fullname="Va_O1_1")
    assert actual == "<i>V</i><sub>O1</sub><sup>1</sup>"


def test_typical_defect_name():
    assert typical_defect_name("Va_O1") is True
    assert typical_defect_name("Va_O") is True

    assert typical_defect_name("split") is False
    assert typical_defect_name("Cu2_split") is False


def test_prettify_names():
    with pytest.raises(ValueError):
        prettify_names({"Va_O": 1, "Va_O1": 1}, None)

    d = {"Va_O1": 1, "Va_O2": 2, "Va_Mg1": 3}

    actual = prettify_names(d, None)
    expected = {'Va_O1': 1, 'Va_O2': 2, "Va_Mg": 3}
    assert actual == expected

    actual = prettify_names(d, "mpl")
    expected = {'$V_{{\\rm Mg}}$': 3, '$V_{{\\rm O}1}$': 1,
                '$V_{{\\rm O}2}$': 2}
    assert actual == expected

    actual = prettify_names(d, "plotly")
    expected = {'<i>V</i><sub>Mg</sub>': 3, '<i>V</i><sub>O1</sub>': 1,
                '<i>V</i><sub>O2</sub>': 2}
    assert actual == expected

    with pytest.raises(ValueError):
        prettify_names(d, "strange_style")
