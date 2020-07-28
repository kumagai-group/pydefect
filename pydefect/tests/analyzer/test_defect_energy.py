# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.analyzer.defect_energy import (
    DefectEnergy, CrossPoints, SingleDefectEnergy, make_defect_energies)


@pytest.fixture
def defect_energy():
    return DefectEnergy(name="Va_O1",
                        charges=[0, 1, 2],
                        energies=[4, 2, -4],
                        corrections=[2, 1, 0])


def test_defect_energy_cross_points(defect_energy):
    actual = defect_energy.cross_points(1, 6)
    expected = CrossPoints([[5.0, 6.0]], [[1.0, -2.0], [6.0, 6.0]])
    assert actual == expected


def test_defect_energy_str(defect_energy):
    expected = """     Va_O1    0       4.0000       2.0000
     Va_O1    1       2.0000       1.0000
     Va_O1    2      -4.0000       0.0000"""
    assert str(defect_energy) == expected


@pytest.fixture
def cross_points():
    inner_cross_points = [[2, 30], [3, 40]]
    boundary_points = [[1, 10], [4, 40]]
    return CrossPoints(inner_cross_points, boundary_points)


def test_cross_points(cross_points):
    assert cross_points.all_sorted_points == [[1, 10], [2, 30], [3, 40], [4, 40]]
    assert cross_points.t_all_sorted_points == [[1, 2, 3, 4], [10, 30, 40, 40]]
    assert cross_points.t_inner_cross_points == [[2, 3], [30, 40]]
    assert cross_points.t_boundary_points == [[1, 4], [10, 40]]
    assert cross_points.charges == [20, 10, 0]
    assert cross_points.charge_list == [(None, 20), (20, 10), (10, 0), (0, None)]


def test_cross_points_str(cross_points):
    expected = """      1.0000      10.0000
      2.0000      30.0000
      3.0000      40.0000
      4.0000      40.0000"""
    assert str(cross_points) == expected


def test_generate_defect_energies():
    defect_energies = [
        SingleDefectEnergy("Va_O1", 0, -1, 1),
        SingleDefectEnergy("Va_Mg1", 0, -2, 2),
        SingleDefectEnergy("Va_O1", 1, 7, 3),
        SingleDefectEnergy("Va_Mg1", -1, -14, 4),
        SingleDefectEnergy("Va_O1", 2, 15, 5),
        SingleDefectEnergy("Va_Mg1", -2, -14, 6)]
    actual = make_defect_energies(defect_energies)
    expected = [DefectEnergy("Va_Mg1", [0, -1, -2], [-2, -14, -14], [2, 4, 6]),
                DefectEnergy("Va_O1", [0, 1, 2], [-1, 7, 15], [1, 3, 5])]
    assert actual == expected

"""
TODO
- Evaluate the crossing points at given Fermi level range.


DONE
- Draw a single defect
- Correction.
"""