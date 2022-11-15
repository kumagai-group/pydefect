# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest
from pydefect.analyzer.defect_energy import CrossPoints
from pydefect.analyzer.transition_levels import make_transition_levels, \
    show_transition_levels, TransitionLevel


@pytest.fixture
def cross_points():
    inner_cross_points = [[2, 30], [3, 40]]
    boundary_points = [[1, 10], [4, 40]]
    return CrossPoints(inner_cross_points, boundary_points)


def test_make_transition_levels(cross_points):
    cross_point_dicts = {"Va_O1": cross_points}
    actual = make_transition_levels(cross_point_dicts)
    expected = [TransitionLevel("Va_O1", [[20, 10], [10, 0]], [30, 40], [2, 3])]
    assert actual == expected


def test_show_transition_levels():
    tl1 = TransitionLevel("Va_O1", [[2, 1], [1, 0]],
                         [1.23456789, 2.23456789], [3.23456789, 4.23456789])
    tl2 = TransitionLevel("Va_Mg1", [[-2, -1], [-1, 0]],
                         [10.23456789, 20.23456789], [30.23456789, 40.23456789])
    print(show_transition_levels([tl1, tl2]))