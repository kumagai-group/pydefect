# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np

from pydefect.util.structure_tools import Distances, Coordination


def test_distances(ortho_conventional):
    distances = Distances(ortho_conventional, center_coord=[0.25, 0.25, 0.25])
    np.testing.assert_almost_equal(distances.distances(), [5.2440442 / 2] * 8)

    distances = Distances(ortho_conventional, center_coord=[0.5, 0.5, 0.5])
    expected = [5.2440442, 3.5, 3.0, 2.5, 3.9051248, 4.3011626, 4.6097722]
    np.testing.assert_almost_equal(distances.distances(), expected)


def test_distances_w_specie(ortho_conventional):
    distances = Distances(ortho_conventional, center_coord=[0.25, 0.26, 0.0])
    expected = [float("inf")] * 4 + [4.030645109656765, 1.9068560512005097,
                                     1.9990247622278219, 3.985737071107425]
    np.testing.assert_equal(distances.distances(remove_self=False, specie="He"),
                            expected)


def test_distances_mapped_atom_idx(ortho_conventional):
    distances = Distances(ortho_conventional, center_coord=[0.0, 0.0, 0.49])
    assert distances.atom_idx_at_center(specie="He") == 4
    distances = Distances(ortho_conventional, center_coord=[0.251, 0.25, 0.25])
    assert distances.atom_idx_at_center(specie="H") is None

    distances = Distances(ortho_conventional, center_coord=[0.0, 0.0, 0.5])
    assert distances.atom_idx_at_center(specie="Li") == None


def test_shortest_distances(ortho_conventional):
    distances = Distances(ortho_conventional, center_coord=[0.5, 0.5, 0.5])
    assert distances.shortest_distance == 2.5


def test_coord_distances(mocker, ortho_conventional):
    mock = mocker.patch("pydefect.util.structure_tools.defaults")
    mock.cutoff_distance_factor = 3.9051248 / 2.5 + 1e-5  # = 1.562
    distances = Distances(ortho_conventional, center_coord=[0.5, 0.5, 0.5])

    actual = distances.coordination
    expected = Coordination({"H": [2.5, 3.0, 3.5], "He": [3.91]},
                            3.905,
                            neighboring_atom_indices={1, 2, 3, 4})
    assert actual == expected
