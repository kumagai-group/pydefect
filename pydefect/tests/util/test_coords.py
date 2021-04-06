# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.util.coords import pretty_coords


def test_pretty_coords():
    assert pretty_coords([0.1, 0.2, 0.3]) == "( 0.100,  0.200,  0.300)"
    assert pretty_coords([1, -2, 0]) == "( 1.000, -2.000,  0.000)"
