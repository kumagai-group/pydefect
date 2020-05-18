# -*- coding: utf-8 -*-

import os
import sys
from math import sqrt
from pathlib import Path

import pytest
from pymatgen import Lattice, IStructure, Structure

sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
# Need the following the share the fixture
from pydefect.tests.input_maker.test_supercell_info import supercell_info


@pytest.fixture(scope="session")
def test_data_files():
    return Path(__file__).parent / "test_data_files"


@pytest.fixture(scope="session")
def simple_cubic():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture(scope="session")
def simple_cubic_2x1x1():
    lattice = Lattice.orthorhombic(2.0, 1.0, 1.0)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H", "H"], coords=coords)


@pytest.fixture(scope="session")
def simple_cubic_2x2x2():
    lattice = Lattice.cubic(2.0)
    coords = [[0.0, 0.0, 0.0],
              [0.5, 0.5, 0.0],
              [0.5, 0.0, 0.5],
              [0.0, 0.5, 0.5],
              [0.0, 0.0, 0.5],
              [0.0, 0.5, 0.0],
              [0.5, 0.0, 0.0],
              [0.5, 0.5, 0.5]]
    return IStructure(lattice=lattice, species=["H"] * 8, coords=coords)


@pytest.fixture(scope="session")
def monoclinic():
    lattice = Lattice.monoclinic(3, 4, 5, 100)
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture(scope="session")
def complex_monoclinic():
    lattice = Lattice.monoclinic(3, 4, 5, 100)
    coords = [[0.0, 0.0, 0.0],
              [0.1, 0.0, 0.0],
              [0.9, 0.0, 0.0],
              [0.2, 0.0, 0.0],
              [0.8, 0.0, 0.0]]
    return IStructure(lattice=lattice,
                          species=["H", "He", "He", "He", "He"],
                          coords=coords)


@pytest.fixture(scope="session")
def rhombohedral():
    lattice = Lattice.rhombohedral(a=1, alpha=45)
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture(scope="session")
def elongated_tetragonal():
    lattice = Lattice.tetragonal(a=1, c=3 * sqrt(2))
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture(scope="session")
def elongated_tetragonal():
    lattice = Lattice.tetragonal(a=1, c=3 * sqrt(2))
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture(scope="session")
def a_centered_orthorhombic():
    lattice = Lattice([[1,  0, 0],
                       [0,  2, 3],
                       [0, -2, 3]])
    coords = [[0.5, 0.8, 0.8],
              [0.0, 0.3, 0.0],
              [0.0, 0.0, 0.3]]

    return IStructure(lattice=lattice, species=["H"] * 3, coords=coords)


@pytest.fixture(scope="session")
def c_centered_monoclinic():
    lattice = Lattice.monoclinic(3, 4, 5, 100)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]]
    return IStructure(lattice=lattice, species=["H", "H"], coords=coords)


@pytest.fixture(scope="session")
def fcc():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]]
    return IStructure(lattice=lattice, species=["H"] * 4, coords=coords)


@pytest.fixture(scope="session")
def bcc():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    return IStructure(lattice=lattice, species=["H"] * 2, coords=coords)


@pytest.fixture
def tetra_close_to_cubic():
    lattice = Lattice.tetragonal(1.001 * 10 / sqrt(2), 10)
    coords = [[0.0, 0.0, 0.0]]
    return IStructure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture
def ortho_conventional():
    lattice = Lattice.orthorhombic(5, 6, 7)
    coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],

        [0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ]
    return IStructure(lattice=lattice, species=["H"] * 4 + ["He"] * 4,
                      coords=coords)


@pytest.fixture
def cubic_supercell():
    lattice = Lattice.cubic(10)
    coords = [
       [0.00, 0.00, 0.00],
       [0.00, 0.00, 0.50],
       [0.00, 0.50, 0.00],
       [0.00, 0.50, 0.50],
       [0.50, 0.00, 0.00],
       [0.50, 0.00, 0.50],
       [0.50, 0.50, 0.00],
       [0.50, 0.50, 0.50],
       [0.00, 0.25, 0.25],
       [0.00, 0.25, 0.75],
       [0.00, 0.75, 0.25],
       [0.00, 0.75, 0.75],
       [0.50, 0.25, 0.25],
       [0.50, 0.25, 0.75],
       [0.50, 0.75, 0.25],
       [0.50, 0.75, 0.75],
       [0.25, 0.00, 0.25],
       [0.25, 0.00, 0.75],
       [0.25, 0.50, 0.25],
       [0.25, 0.50, 0.75],
       [0.75, 0.00, 0.25],
       [0.75, 0.00, 0.75],
       [0.75, 0.50, 0.25],
       [0.75, 0.50, 0.75],
       [0.25, 0.25, 0.00],
       [0.25, 0.25, 0.50],
       [0.25, 0.75, 0.00],
       [0.25, 0.75, 0.50],
       [0.75, 0.25, 0.00],
       [0.75, 0.25, 0.50],
       [0.75, 0.75, 0.00],
       [0.75, 0.75, 0.50],
       [0.25, 0.00, 0.00],
       [0.25, 0.00, 0.50],
       [0.25, 0.50, 0.00],
       [0.25, 0.50, 0.50],
       [0.75, 0.00, 0.00],
       [0.75, 0.00, 0.50],
       [0.75, 0.50, 0.00],
       [0.75, 0.50, 0.50],
       [0.25, 0.25, 0.25],
       [0.25, 0.25, 0.75],
       [0.25, 0.75, 0.25],
       [0.25, 0.75, 0.75],
       [0.75, 0.25, 0.25],
       [0.75, 0.25, 0.75],
       [0.75, 0.75, 0.25],
       [0.75, 0.75, 0.75],
       [0.00, 0.00, 0.25],
       [0.00, 0.00, 0.75],
       [0.00, 0.50, 0.25],
       [0.00, 0.50, 0.75],
       [0.50, 0.00, 0.25],
       [0.50, 0.00, 0.75],
       [0.50, 0.50, 0.25],
       [0.50, 0.50, 0.75],
       [0.00, 0.25, 0.00],
       [0.00, 0.25, 0.50],
       [0.00, 0.75, 0.00],
       [0.00, 0.75, 0.50],
       [0.50, 0.25, 0.00],
       [0.50, 0.25, 0.50],
       [0.50, 0.75, 0.00],
       [0.50, 0.75, 0.50],
    ]
    return Structure(lattice=lattice, species=["H"] * 32 + ["He"] * 32,
                     coords=coords)
