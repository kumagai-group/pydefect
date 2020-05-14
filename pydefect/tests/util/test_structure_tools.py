# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
import numpy as np

from pydefect.util.structure_tools import Distances
from pymatgen import Lattice, Structure


def test_distances(monoclinic_supercell):
    distances = Distances(monoclinic_supercell, coord=[0.25, 0.25, 0.25])
    np.testing.assert_almost_equal(distances.distances(), [5.2440442 / 2] * 8)

    distances = Distances(monoclinic_supercell, coord=[0.5, 0.5, 0.5])
    expected = [5.2440442, 3.5, 3.0, 2.5, 3.9051248, 4.3011626, 4.6097722]
    np.testing.assert_almost_equal(distances.distances(), expected)


def test_shortest_distances(monoclinic_supercell):
    distances = Distances(monoclinic_supercell, coord=[0.5, 0.5, 0.5])
    assert distances.shortest_distance == 2.5


def test_coord_distances(monoclinic_supercell):
    distances = Distances(monoclinic_supercell, coord=[0.5, 0.5, 0.5])
    cutoff_distance_factor = 3.9051248 / 2.5 + 1e-5

    actual = distances.coord_distances(cutoff_distance_factor)
    expected = {"H": [2.5, 3.0, 3.5], "He": [3.91]}
    assert actual == expected


"""
TODO

DONE
- Return distances to atoms in the given structure from the given coord.
- Return the shortest distance from the given coord to atoms
- Return the coordination distances from the given coord to atoms
  The multiplication_factor to the shortest distance is given.
"""