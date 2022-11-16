# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest
from pydefect.analyzer.defect_energy import CrossPoints
from pydefect.analyzer.transition_levels import make_transition_levels, \
    TransitionLevel, TransitionLevels
from vise.tests.helpers.assertion import assert_msonable


@pytest.fixture
def cross_points():
    inner_cross_points = [[2, 30], [3, 40]]
    boundary_points = [[1, 10], [4, 40]]
    return CrossPoints(inner_cross_points, boundary_points)


@pytest.fixture
def transition_levels():
    return TransitionLevels(
        [TransitionLevel("Va_O1", [[20, 10], [10, 0]], [30, 40], [2, 3])],
        3.0, -0.1, 3.0)


def test_transition_levels_msonable(transition_levels):
    assert_msonable(transition_levels)


def test_make_transition_levels(cross_points, transition_levels):
    cross_point_dicts = {"Va_O1": cross_points}
    actual = make_transition_levels(cross_point_dicts, 3.0, -0.1, 3.0)
    assert actual == transition_levels


def test_show_transition_levels():
    tl1 = TransitionLevel("Va_O1", [[2, 1], [1, 0]],
                          [1.23456789, 2.23456789], [3.23456789, 4.23456789])
    tl2 = TransitionLevel("Va_Mg1", [[-2, -1], [-1, 0]],
                          [10.2345678, 20.2345678], [30.2345678, 40.2345678])
    tls = TransitionLevels([tl1, tl2],  3.0, -0.1, 3.0)
    actual = tls.__str__()
    expected = """vbm: 0.00, cbm: 3.00, supercell vbm: -0.10, supercell_cbm: 3.00
name    charges      Fermi level    Formation energy
------  ---------  -------------  ------------------
Va_O1   2 | 1              3.235               1.235
        1 | 0              4.235               2.235
Va_Mg1  -2 | -1           30.235              10.235
        -1 | 0            40.235              20.235"""

    assert actual == expected