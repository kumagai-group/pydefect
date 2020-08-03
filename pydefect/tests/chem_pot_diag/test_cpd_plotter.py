# -*- coding: utf-8 -*-

import numpy as np
import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo
from pydefect.chem_pot_diag.cpd_plotter import (
    ChemPotDiagMpl2DMplPlotter, transpose, sort_coords,
    ChemPotDiagMpl3DMplPlotter, ChemPotDiagPlotly2DMplPlotter)
from pymatgen.core.composition import Composition


@pytest.fixture
def cpd_plot_info_2d():
    energies = {Composition("H"): 0.0,
                Composition("O"): 1.0,
                Composition("H4O2"): -4.0}
    cpd = ChemPotDiag(energies, target=Composition("H2O"))
    return CpdPlotInfo(cpd, min_range=-10)


@pytest.mark.skipif(False, reason="")
def test_cpd_2d_draw(cpd_plot_info_2d):
    plotter = ChemPotDiagMpl2DMplPlotter(cpd_plot_info_2d)
    plotter.draw_diagram().show()


@pytest.fixture
def cpd_3d_info():
    energies = {Composition("H"): 0.0,
                Composition("O"): 1.0,
                Composition("H2O"): -2.0,
                Composition("Mg"): 0.0,
                Composition("MgO"): -3.0}
    cpd = ChemPotDiag(energies, target=Composition("H2O"))
    return CpdPlotInfo(cpd)


@pytest.mark.skipif(False, reason="")
def test_cpd_3d_draw(cpd_3d_info):
    plotter = ChemPotDiagMpl3DMplPlotter(cpd_3d_info)
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


def test_plotly_2d(cpd_plot_info_2d):
    plotter = ChemPotDiagPlotly2DMplPlotter(cpd_plot_info_2d)
    fig = plotter.figure()
    fig.show()
