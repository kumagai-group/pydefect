# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.cli.vasp.utils.main_vasp_util_functions import make_parchg_dir
from vise.input_set.incar import ViseIncar


def test_make_parchg_dir(tmpdir, mocker):
    print(tmpdir)
    tmpdir.chdir()
    Path("WAVECAR").touch()
    Path("INCAR").write_text("ALGO  =  Normal")
    Path("POSCAR").write_text("a")
    Path("POTCAR").write_text("b")
    Path("KPOINTS").write_text("c")

    band_edge_states = mocker.Mock(spec=BandEdgeStates, autospec=True)
    band_edge_states.band_indices_from_vbm_to_cbm = [1, 2]

    with pytest.raises(FileNotFoundError):
        make_parchg_dir(band_edge_states=band_edge_states)

    Path("WAVECAR").write_text("d")
    make_parchg_dir(band_edge_states=band_edge_states)

    assert Path("parchg/POSCAR").read_text() == "a"
    assert Path("parchg/POTCAR").read_text() == "b"
    assert Path("parchg/KPOINTS").read_text() == "c"
    assert Path("parchg/WAVECAR").read_text() == "d"

    actual = ViseIncar.from_file("parchg/INCAR")
    expected = {"ALGO": "Normal", "LPARD": True, "LSEPB": True, "KPAR": 1,
                "IBAND": "2 3"}
    assert actual == expected
