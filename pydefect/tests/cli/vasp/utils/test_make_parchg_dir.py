# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeStates, BandEdgeState
from pydefect.cli.vasp.utils.make_parchg_dir import make_parchg_dir
from pymatgen import Lattice, Structure


def test_make_parchg_dir(tmpdir, mocker):
    print(tmpdir)
    tmpdir.chdir()
    Path("WAVECAR").touch()
    band_edge_states = mocker.Mock(spec=BandEdgeStates, autospec=True)
    band_edge_states.band_indices_from_vbm_to_cbm = [1, 2]

    with pytest.raises(FileNotFoundError):
        make_parchg_dir(band_edge_states=band_edge_states)

    s = Structure(Lattice.cubic(1.0), species=["H"], coords=[[0.0]*3])
    s.to(filename="CONTCAR", fmt="POSCAR")
    Path("WAVECAR").write_text("a")
    make_parchg_dir(band_edge_states=band_edge_states)

