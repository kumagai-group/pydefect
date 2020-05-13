# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest

from pymatgen import Lattice, Structure

from pydefect.input_maker.supercell import Supercell, Supercells


@pytest.fixture
def simple_cubic_supercell():
    lattice = Lattice.orthorhombic(2.0, 1.0, 1.0)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
    results = Structure(lattice=lattice, species=["H", "H"], coords=coords)

    return results


def test_supercell(simple_cubic, simple_cubic_supercell):
    matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    supercell = Supercell(input_structure=simple_cubic, matrix=matrix)
    assert supercell.structure == simple_cubic_supercell

    average = (2 + 1 + 1) / 3
    expected = (abs(2 - average) + abs(1 - average) * 2) / 3 / average
    assert supercell.isotropy == expected


def test_supercell_average_angle(monoclinic):
    matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    supercell = Supercell(input_structure=monoclinic, matrix=matrix)
    assert supercell.average_angle == (90 + 90 + 100) / 3


def test_supercells(simple_cubic):
    supercells = Supercells(input_structure=simple_cubic,
                            min_num_atoms=64,
                            max_num_atoms=100)
    expected = simple_cubic * [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
    assert isinstance(supercells.isotropic_supercell, Supercell)
    assert supercells.isotropic_supercell.structure == expected

    # supercells = Supercells(input_structure=simple_cubic,
    #                         min_num_atoms=65,
    #                         max_num_atoms=100)
    # expected = simple_cubic * [[5, 0, 0], [0, 5, 0], [0, 0, 5]]

    # assert supercells.isotropic_supercell == expected


"""
TODO
- 

DONE
- Input of structure and trans_mat returns supercell
- Return isotropy
- Return average angle
- Recommend the least isotropic supercell of simple cubic within a given number of atom range.
  When the isotropy is the same, smallest supercell is returned.
- Check max_num_atoms

REJECT
+ Allow to write the supercell down a file.
    -> Structure already has to, which should be used.
"""