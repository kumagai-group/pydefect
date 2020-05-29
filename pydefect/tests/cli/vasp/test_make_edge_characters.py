# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest
from pymatgen import Spin, Structure, Lattice
from pymatgen.io.vasp import Vasprun, Outcar, Procar

from pydefect.cli.vasp.make_edge_characters import MakeEdgeCharacters, \
    calc_participation_ratio, calc_orbital_character


@pytest.fixture
def make_edge_characters(mocker):
    mock_procar = mocker.Mock(spec=Procar, autospec=True)
    mock_vasprun = mocker.Mock(spec=Vasprun, autospec=True)
    mock_outcar = mocker.Mock(spec=Outcar, autospec=True)
#    mock_ratio = mocker.patch("pydefect.cli.vasp.make_edge_characters.calc_participation_ratio")
    mock_procar.data = {Spin.up: np.array(
        # 1st atom
        # first k-index
        [[[[0.0] * 15 + [0.1]] * 3,   # 1st band
          [[0.1] * 16, [0.2] * 16, [1.0] * 16],
          [[0.0] * 15 + [0.1]] * 3],  # 2nd band
         # 2nd k-index      2nd k-index
         [[[0.0] * 15 + [0.1]] * 3,
          [[0.0] * 15 + [0.1]] * 3,
          [[0.0] * 15 + [0.1]] * 3]]),
    Spin.down: np.array(
        # 1st atom
        # first k-index
        [[[[0.0] * 15 + [0.1]] * 3,   # 1st band
          [[0.1] * 16, [0.2] * 16, [1.0] * 16],
          [[0.0] * 15 + [0.1]] * 3],  # 2nd band
         # 2nd k-index      2nd k-index
         [[[0.0] * 15 + [0.1]] * 3,
          [[0.0] * 15 + [0.1]] * 3,
          [[0.0] * 15 + [0.1]] * 3]])}
    mock_vasprun.final_structure = Structure(Lattice.cubic(1), species=["H", "H", "He"], coords=[[0] * 3] * 3)
    mock_vasprun.eigenvalues = {Spin.up: np.array([[[-3.01, 1.],
                                                    [-2.90, 1.],  # vbm
                                                    [ 8.01, 0.]],
                                                   [[-3.01, 1.],
                                                    [-2.99, 1.],
                                                    [ 7.90, 0.]]]), # cbm
                    # shift one band for testing magnetization
                                Spin.down: np.array([[[-2.99, 1.],
                                                      [ 7.99, 0.],
                                                      [10.00, 0.]],
                                                     [[-3.01, 1.],
                                                      [ 7.99, 0.],
                                                      [10.00, 0.]]])}
    mock_outcar.nelect = 3.0
    mock_outcar.total_mag = 1

    return MakeEdgeCharacters(mock_procar, mock_vasprun, mock_outcar, [0, 1])


def test_edge_characters(make_edge_characters):
    edge_characters = make_edge_characters.edge_characters
    assert edge_characters[0].vbm == -2.90
    assert edge_characters[0].hob_bottom_e == -2.99
    assert edge_characters[0].cbm == 7.90
    assert edge_characters[0].lub_top_e == 8.01

    assert edge_characters[1].vbm == -2.99
    assert edge_characters[1].hob_bottom_e == -3.01
    assert edge_characters[1].cbm == 7.99
    assert edge_characters[1].lub_top_e == 7.99


def test_calc_participation_ratio():
    # only a single atom
    actual = calc_participation_ratio(orbitals={Spin.up: np.array(
        # first k-index
        [[[[0.0, 0.0, 0.1], [0.0, 0.0, 0.1], [0.0, 0.0, 0.1]],   # 1st band
          [[0.2, 0.2, 0.2], [0.2, 0.2, 0.2], [0.3, 0.3, 0.3]]],  # 2nd band
         # 2nd k-index
         [[[0.0, 0.0, 0.1], [0.0, 0.0, 0.1], [0.0, 0.0, 0.1]],
          [[0.0, 0.0, 0.1], [0.0, 0.0, 0.1], [0.0, 0.0, 0.1]]]])},
        spin=Spin.up, kpt_index=0, band_index=1, atom_indices=[0, 1])
    np.testing.assert_almost_equal(actual, 1.2 / (1.2 + 0.9))


def test_calc_orbital_character_azimuthal_false():
    spd_orbitals = {Spin.up: np.array(
        # 1st atom
        # first k-index
        [[[[0.0, 0.0, 0.1], [0.0, 0.0, 0.1], [0.0, 0.0, 0.1]],   # 1st band
          [[0.1, 0.1, 0.1], [0.2, 0.3, 0.4], [1.0, 1.1, 1.2]]],  # 2nd band
         # 2nd k-index      2nd k-index
         [[[0.0, 0.0, 0.1], [0.0, 0.0, 0.1], [0.0, 0.0, 0.1]],
          [[0.0, 0.0, 0.1], [0.0, 0.0, 0.1], [0.0, 0.0, 0.1]]]])}

    actual = calc_orbital_character(
        orbitals=spd_orbitals,
        structure=Structure(Lattice.cubic(1),
                            species=["H", "H", "He"], coords=[[0] * 3] * 3),
        spin=Spin.up, kpt_index=0, band_index=1)

    assert actual == {"H": [0.3, 0.4, 0.5], "He": [1.0, 1.1, 1.2]}


def test_calc_orbital_character_azimuthal_true():
    azimuthal_orbitals = {Spin.up: np.array(
        # 1st atom
        # first k-index
        [[[[0.0] * 15 + [0.1]] * 3,   # 1st band
          [[0.1] * 16, [0.2] * 16, [1.0] * 16]],  # 2nd band
         # 2nd k-index      2nd k-index
         [[[0.0] * 15 + [0.1]] * 3,
          [[0.0] * 15 + [0.1]] * 3]])}

    actual = calc_orbital_character(
        orbitals=azimuthal_orbitals,
        structure=Structure(Lattice.cubic(1),
                            species=["H", "H", "He"], coords=[[0] * 3] * 3),
        spin=Spin.up, kpt_index=0, band_index=1)

    assert actual == {"H": [0.3, 0.9, 1.5, 2.1], "He": [1.0, 3.0, 5.0, 7.0]}


"""
TODO
-

DONE
"""