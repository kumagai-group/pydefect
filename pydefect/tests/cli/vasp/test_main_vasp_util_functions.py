# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.analyzer.calc_results import CalcResults
from pydefect.cli.vasp.main_vasp_util_functions import make_parchg_dir, \
    make_refine_defect_poscar, calc_charge_state, \
    calc_grids, make_defect_charge_info_main, make_defect_entry_main
from pymatgen.core import Structure
from vise.input_set.incar import ViseIncar
import numpy as np

_filepath = "pydefect.cli.vasp.main_vasp_util_functions"


def test_calc_charge_state(mocker):
    mock_poscar = mocker.patch(f"{_filepath}.Poscar")
    mock_potcar = mocker.patch(f"{_filepath}.Potcar")
    mock_incar = mocker.patch(f"{_filepath}.Incar")
    mock_get_charge_state = mocker.patch(f"{_filepath}.get_defect_charge_state")
    mock_get_charge_state.return_value = 0

    args = Namespace(dir=Path("Va_O1_0"))
    calc_charge_state(args)

    mock_poscar.from_file.assert_called_once_with(Path("Va_O1_0/POSCAR"))
    mock_potcar.from_file.assert_called_once_with(Path("Va_O1_0/POTCAR"))
    mock_incar.from_file.assert_called_once_with(Path("Va_O1_0/INCAR"))
    mock_get_charge_state.assert_called_once_with(
        mock_poscar.from_file.return_value,
        mock_potcar.from_file.return_value,
        mock_incar.from_file.return_value)


def test_make_defect_entry_main(mocker):
    mock_charge_state = mocker.patch(f"{_filepath}.calc_charge_state")
    mock_charge_state.return_value = 0
    mock_structure = mocker.patch(f"{_filepath}.Structure")
    mock_make_defect_entry = mocker.patch(f"{_filepath}.make_defect_entry")
    mock_perfect = mocker.Mock()

    args = Namespace(dir=Path("Va_O1_0"), name="Va_O1", perfect=mock_perfect)
    make_defect_entry_main(args)
    mock_charge_state.assert_called_once_with(args)
    mock_structure.from_file.assert_called_once_with(Path("Va_O1_0/POSCAR"))
    mock_make_defect_entry.assert_called_once_with(
        name="Va_O1",
        charge=0,
        perfect_structure=mock_perfect,
        defect_structure=mock_structure.from_file.return_value)


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

    mock_loadfn = mocker.patch(f"{_filepath}.loadfn", side_effect=side_effect)
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


def test_calc_grids(mocker):
    mock_grids = mocker.patch(f"{_filepath}.Grids")
    args = Namespace(chgcar="CHGCAR")
    calc_grids(args)
    mock_grids.from_chgcar.assert_called_once_with("CHGCAR")
    mock_grids.from_chgcar.return_value.dump.assert_called_once_with()


def test_make_defect_charge_info_main(mocker):
    mock_chgcar = mocker.patch(f"{_filepath}.Chgcar")
    m_chgcar = mock_chgcar.from_file.return_value
    mock_make_charge_info = mocker.patch(f"{_filepath}.make_defect_charge_info")
    mock_charge_info = mock_make_charge_info.return_value
    mock_grids = mocker.Mock()
    args = Namespace(parchgs=["PARCHG.0189.ALLK"],
                     grids=mock_grids,
                     bin_interval=0.1)
    make_defect_charge_info_main(args)

    mock_make_charge_info.assert_called_once_with([m_chgcar], [188], 0.1, mock_grids)
    mock_chgcar.from_file.assert_called_once_with("PARCHG.0189.ALLK")
    mock_charge_info.to_json_file.assert_called_once_with()
