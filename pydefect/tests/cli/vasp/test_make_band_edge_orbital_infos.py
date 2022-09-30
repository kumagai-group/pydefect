# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import numpy as np
from pydefect.analyzer.band_edge_states import OrbitalInfo, \
    BandEdgeOrbitalInfos
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.cli.vasp.make_band_edge_orbital_infos import \
    make_band_edge_orbital_infos
from pymatgen.core import Structure, Lattice
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Procar, Vasprun
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


def test_make_band_edge_orbital_infos(mocker):
    mock_procar = mocker.Mock(spec=Procar, autospec=True)
    mock_vasprun = mocker.Mock(spec=Vasprun, autospec=True)
    mock_str_info = mocker.Mock(spec=DefectStructureInfo, autospec=True)

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
    mock_vasprun.efermi = 20.0
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

    mock_str_info.neighbor_atom_indices = [0]
    actual = make_band_edge_orbital_infos(
        mock_procar, mock_vasprun, vbm=0.0, cbm=5.0, str_info=mock_str_info)

    expected = BandEdgeOrbitalInfos(
        orbital_infos=[[
            [OrbitalInfo(energy=-2.9, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=1.0, participation_ratio=1.0),
             OrbitalInfo(energy=8.01, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0, participation_ratio=1.0)],
            [OrbitalInfo(energy=-2.99, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=1.0, participation_ratio=1.0),
             OrbitalInfo(energy=7.90, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0, participation_ratio=1.0)]],
           [[OrbitalInfo(energy=7.99, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0, participation_ratio=1.0),
             OrbitalInfo(energy=10.00, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0, participation_ratio=1.0)],
            [OrbitalInfo(energy=8.00, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0, participation_ratio=1.0),
             OrbitalInfo(energy=10.00, orbitals={"H": [1.0, 0.0, 0.0, 0.0], "He": [0.0, 0.0, 0.0, 0.0]}, occupation=0.0, participation_ratio=1.0)]]],
        kpt_coords=[(0.0, 0.0, 0.0)], kpt_weights=[1.0], lowest_band_index=1,
    fermi_level=20.0)
    assert_dataclass_almost_equal(actual, expected, check_is_subclass=True)

