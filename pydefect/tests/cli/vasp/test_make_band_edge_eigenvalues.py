# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pydefect.cli.vasp.make_band_edge_eigenvalues import \
    make_band_edge_eigenvalues
from pymatgen import Spin
from pymatgen.io.vasp import Vasprun


@pytest.fixture
def band_edge_eigenvalues(mocker):
    mock_vasprun = mocker.Mock(spec=Vasprun, autospec=True)
    mock_vasprun.actual_kpoints = [[0.0, 0.0, 0.0]]
    mock_defaults = mocker.patch("pydefect.cli.vasp.make_band_edge_eigenvalues.defaults")
    mock_defaults.eigval_range = 3.0
    mock_vasprun.eigenvalues = {Spin.up: np.array([[[-3.01, 1.],
                                           [-2.99, 1.],
                                           [ 8.01, 0.],
                                           [10.00, 0.]],
                                          [[-3.01, 1.],
                                           [-2.99, 1.],
                                           [ 7.99, 0.],
                                           [10.00, 0.]]]),
                                Spin.down: np.array([[[-3.01, 1.],
                                           [-2.99, 1.],
                                           [ 7.99, 0.],  # only change here from Spin.up
                                           [10.00, 0.]],
                                          [[-3.01, 1.],
                                           [-2.99, 1.],
                                           [ 7.99, 0.],
                                           [10.00, 0.]]])}
    return make_band_edge_eigenvalues(mock_vasprun, vbm=0.0, cbm=5.0)


def test_make_band_edge_eigenvalues_kpt_coords(band_edge_eigenvalues):
    assert band_edge_eigenvalues.kpt_coords == [(0.0, 0.0, 0.0)]


def test_make_band_edge_eigenvalues_eigenvalues(band_edge_eigenvalues):
    actual = band_edge_eigenvalues.energies_and_occupations
    expected = [[[[-2.99, 1.], [8.01, 0.]], [[-2.99, 1.], [7.99, 0.]]],
                [[[-2.99, 1.], [7.99, 0.]], [[-2.99, 1.], [7.99, 0.]]]]
    np.testing.assert_array_almost_equal(actual, expected)


def test_make_band_edge_eigenvalues_lowest_band_index(band_edge_eigenvalues):
    assert band_edge_eigenvalues.lowest_band_index == 1
