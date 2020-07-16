# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues
from pydefect.analyzer.eigenvalue_plotter import EigenvalueMplPlotter


@pytest.fixture
def eigenvalue_plotter():
    eig = BandEdgeEigenvalues(energies_and_occupations=[
        [[[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]], [[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]]],
        [[[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]], [[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]]],
                                                        ],
        kpt_coords=[(0.0, 0.0, 0.0), (0.25, 0.0, 0.0)],
        lowest_band_index=10)
    return EigenvalueMplPlotter(title="test", band_edge_eigenvalues=eig,
                                supercell_vbm=0.1, supercell_cbm=0.9)


def test_plot(eigenvalue_plotter):
    eigenvalue_plotter.construct_plot()
    eigenvalue_plotter.plt.show()


def test_plot_with_actual_file(eigenvalue_plotter):
    eig = loadfn(Path(__file__).parent / "band_edge_eigenvalues.json")
    plotter = EigenvalueMplPlotter(title="test", band_edge_eigenvalues=eig,
                                   supercell_vbm=0.5, supercell_cbm=2.6)
    plotter.construct_plot()
    plotter.plt.show()

