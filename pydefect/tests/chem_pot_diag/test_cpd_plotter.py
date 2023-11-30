# -*- coding: utf-8 -*-

import numpy as np
import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, \
    TargetVertex
from pydefect.chem_pot_diag.cpd_plotter import (
    ChemPotDiag2DMplPlotter, transpose, sort_coords,
    ChemPotDiag3DMplPlotter)


@pytest.fixture
def cpd_2d():
    polygons = {"H": [[0.0, -3.3], [0.0, -3.0]],
                "O": [[-3.3, 0.0], [-3.0, 0.0]],
                "HO": [[-3.0, 0.0], [0.0, -3.0]]}
    return ChemPotDiag(
        vertex_elements=["H", "O"],
        polygons=polygons,
        target="HO",
        target_vertices_dict=
        {"A": TargetVertex({"H": 0.0, "O": -3.0}, ["H"], None),
         "B": TargetVertex({"H": -3.0, "O": 0.0}, ["O"], None)})


@pytest.mark.skipif(False, reason="")
def test_cpd_2d_draw(cpd_2d):
    plotter = ChemPotDiag2DMplPlotter(cpd_2d)
    plotter.draw_diagram().show()


@pytest.fixture
def cpd_3d():
    polygons = {"H": [[0.0, 0.0, -3.3], [0.0, 0.0, -3.0], [0.0, -3.0, 0.0], [0.0, -3.3, 0.0], [0.0, -3.3, -3.3]],
                "O": [[0.0, 0.0, -3.3], [0.0, 0.0, -3.0], [-3.0, 0.0, 0.0], [-3.3, 0.0, 0.0], [-3.3, 0.0, -3.3]],
                "N": [[-3.3, 0.0, 0.0], [-3.0, 0.0, 0.0], [0.0, -3.0, 0.0], [0.0, -3.3, 0.0], [-3.3, -3.3, 0.0]],
                "HON":  [[0.0, 0.0, -3.0], [0.0, -3.0, 0.0], [-3.0, 0.0, 0.0]]}
    return ChemPotDiag(
        vertex_elements=["H", "O", "N"],
        polygons=polygons,
        target="HON",
        target_vertices_dict={
            "A": TargetVertex({"H": 0.0, "O": 0.0, "N": -3.0}, ["H", "O"], None),
            "B": TargetVertex({"H": 0.0, "O": -3.0, "N": 0.0}, ["H", "N"], None),
            "C": TargetVertex({"H": -3.0, "O": 0.0, "N": 0.0}, ["O", "N"], None)})


@pytest.mark.skipif(False, reason="")
def test_cpd_3d_draw(cpd_3d):
    plotter = ChemPotDiag3DMplPlotter(cpd_3d)
    plotter.draw_diagram().show()


def test_transpose():
    target_list = [[1, 2], [3, 4]]
    actual = transpose(target_list)
    expected = [[1, 3], [2, 4]]
    assert actual == expected


def test_sort_coords():
    coords = np.array([[3, 2, -1], [-1, 2, 0], [-6, -1, 4], [1, -3, 3]])
    expected = np.array([[1, -3, 3], [-6, -1, 4], [-1, 2, 0], [3, 2, -1]])
    np.testing.assert_array_equal(sort_coords(coords), expected)
