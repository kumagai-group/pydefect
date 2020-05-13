# -*- coding: utf-8 -*-

from pathlib import Path

import pytest
from pymatgen import Lattice, Structure


@pytest.fixture(scope="session")
def test_data_files():
    return Path(__file__).parent / "test_data_files"


@pytest.fixture(scope="session")
def simple_cubic():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    results = Structure(lattice=lattice, species=["H"], coords=coords)

    return results


@pytest.fixture(scope="session")
def monoclinic():
    lattice = Lattice.monoclinic(3, 4, 5, 100)
    coords = [[0.0, 0.0, 0.0]]
    structure = Structure(lattice=lattice, species=["H"], coords=coords)
    return structure


@pytest.fixture(scope="session")
def complex_monoclinic():
    lattice = Lattice.monoclinic(3, 4, 5, 100)
    coords = [[0.0, 0.0, 0.0],
              [0.1, 0.0, 0.0],
              [0.9, 0.0, 0.0],
              [0.2, 0.0, 0.0],
              [0.8, 0.0, 0.0]]
    structure = Structure(lattice=lattice,
                          species=["H", "He", "He", "He", "He"],
                          coords=coords)
    return structure

