# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.analyzer.band_edge_states import OrbitalInfo, \
    BandEdgeOrbitalInfos
from pydefect.cli.vasp.make_band_edge_orbital_characters import \
    make_band_edge_orbital_characters
from pymatgen import Spin, Structure, Lattice
from pymatgen.io.vasp import Procar, Vasprun, Outcar
from pytest_mock import mocker
import numpy as np
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


def test_make_band_edge_orbital_characters(mocker):
    mock_procar = mocker.Mock(spec=Procar, autospec=True)
    mock_vasprun = mocker.Mock(spec=Vasprun, autospec=True)
    mock_defaults = mocker.patch("pydefect.cli.vasp.make_band_edge_eigenvalues.defaults")
    mock_defaults.eigval_range = 3.0

    mock_vasprun.actual_kpoints = [[0.0, 0.0, 0.0]]
    mock_vasprun.actual_kpoints_weights = [1.0]
    mock_vasprun.final_structure = Structure(
        Lattice.cubic(1), species=["H", "H", "He"], coords=[[0] * 3] * 3)
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
                                                      [ 8.00, 0.],
                                                      [10.00, 0.]]])}
    # s: 1, p: 3, d: 5, f: 7 = 16 orbitals
    mock_procar.data = {Spin.up: np.array(
       [[[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1st kpt, 0th band, 1st atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1st kpt, 0th band, 2nd atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]], # 1st kpt, 0th band, 3rd atom
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1st kpt, 1st band, 1st atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]],
        [[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], # 2nd kpt, 0th band, 1st atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]]),
        Spin.down: np.array(
       [[[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1st kpt, 0th band, 1st atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1st kpt, 0th band, 2nd atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]], # 1st kpt, 0th band, 3rd atom
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # 1st kpt, 1st band, 1st atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]],
        [[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], # 2nd kpt, 0th band, 1st atom
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]],
         [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]]]),
    }

    actual = make_band_edge_orbital_characters(
        mock_procar, mock_vasprun, vbm=0.0, cbm=5.0)

    expected = BandEdgeOrbitalInfos(
        orbital_infos=[[
            [OrbitalInfo(energy=-2.9, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=1.0),
             OrbitalInfo(energy=8.01, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0)],
            [OrbitalInfo(energy=-2.99, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=1.0),
             OrbitalInfo(energy=7.90, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0)]],
           [[OrbitalInfo(energy=7.99, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0),
             OrbitalInfo(energy=10.00, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0)],
            [OrbitalInfo(energy=8.00, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0),
             OrbitalInfo(energy=10.00, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0)]]],
        kpt_coords=[(0.0, 0.0, 0.0)], kpt_weights=[1.0], lowest_band_index=1)
    assert_dataclass_almost_equal(actual, expected)

