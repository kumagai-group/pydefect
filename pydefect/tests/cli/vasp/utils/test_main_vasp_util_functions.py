# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.analyzer.calc_results import CalcResults
from pydefect.cli.vasp.utils.main_vasp_util_functions import make_parchg_dir, \
    make_refine_defect_poscar
from pymatgen import Structure
from vise.input_set.incar import ViseIncar
import numpy as np


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

    calc_results = mocker.Mock(spec=CalcResults, autospec=True)

    def side_effect(key):
        if str(key) == "band_edge_states.json":
            return band_edge_states
        elif str(key) == "calc_results.json":
            return calc_results
        else:
            raise ValueError

    mock_loadfn = mocker.patch(
        "pydefect.cli.vasp.utils.main_vasp_util_functions.loadfn",
        side_effect=side_effect)
    args = Namespace(dir=Path(tmpdir))

    with pytest.raises(FileNotFoundError):
        make_parchg_dir(args)

    Path("WAVECAR").write_text("d")
    make_parchg_dir(args)

    assert Path("parchg/POSCAR").read_text() == "a"
    assert Path("parchg/POTCAR").read_text() == "b"
    assert Path("parchg/KPOINTS").read_text() == "c"
    assert Path("parchg/WAVECAR").read_text() == "d"

    actual = ViseIncar.from_file("parchg/INCAR")
    expected = {"ALGO": "Normal", "LPARD": True, "LSEPB": True, "KPAR": 1,
                "IBAND": "2 3"}
    assert actual == expected


def test_refine_poscar(tmpdir, mocker, before_refine, after_refine):
    print(tmpdir)
    tmpdir.chdir()
    mock_defect_entry = mocker.Mock()
    mock_defect_entry.anchor_atom_index = 1
    mock_defect_entry.anchor_atom_coords = np.array([0.0, 0.5, 0.5])

    args = Namespace(structure=before_refine,
                     defect_entry=mock_defect_entry,
                     poscar_name="refined_POSCAR")
    make_refine_defect_poscar(args)

    assert Structure.from_file("refined_POSCAR") == after_refine

