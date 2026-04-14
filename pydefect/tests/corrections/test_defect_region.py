# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.
import numpy as np
import pytest

from pydefect.corrections.defect_region import calc_max_sphere_radius, \
    FixedDistanceDefectRegion, HalfMaxFaceDistanceDefectRegion


def test_calc_max_sphere_radius_simple_cubic():
    lattice = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    assert calc_max_sphere_radius(lattice) == pytest.approx(0.5)

def test_calc_max_sphere_radius_rectangular():
    lattice = np.array([
        [2.0, 0.0, 0.0],
        [0.0, 3.0, 0.0],
        [0.0, 0.0, 4.0],
    ])
    assert calc_max_sphere_radius(lattice) == pytest.approx(2.0)

def test_calc_max_sphere_radius_non_orthogonal():
    sin60 = np.sqrt(3) / 2.0
    lattice = np.array([
        [1.0, 0.0, 0.0],
        [0.5, sin60, 0.0],
        [0.0, 0.0, 1.0],
    ])
    assert calc_max_sphere_radius(lattice) == pytest.approx(0.5)


def test_fixed_distance_defect_region_returns_radius():
    r = 1.23
    region = FixedDistanceDefectRegion(r)
    assert region.defect_region_radius() == pytest.approx(r)

    assert region.defect_region_radius(np.eye(3)) == pytest.approx(r)


def test_half_min_face_distance_defect_region_cube_lattice():
    lattice = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
    ratio = 0.5
    region = HalfMaxFaceDistanceDefectRegion(ratio)

    expected = calc_max_sphere_radius(lattice) * ratio
    assert region.defect_region_radius(lattice) == pytest.approx(expected)