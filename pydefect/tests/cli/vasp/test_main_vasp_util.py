# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.cli.vasp.main_vasp_util import parse_args_main_vasp_util


def test_calc_charge_state():
    parsed_args = parse_args_main_vasp_util(
        ["ccs", "-d", "Va_O1_0"])
    expected = Namespace(
        dir=Path("Va_O1_0"),
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_defect_entry(mocker):
    mock_structure = mocker.patch("pydefect.cli.vasp.main_vasp_util.Structure")
    parsed_args = parse_args_main_vasp_util(
        ["de", "-d", "Va_O1_0", "-n", "Va_O1", "-p", "POSCAR"])
    expected = Namespace(
        dir=Path("Va_O1_0"),
        name="Va_O1",
        perfect=mock_structure.from_file.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_structure.from_file.assert_called_once_with("POSCAR")


def test_make_parchg_dir():
    parsed_args = parse_args_main_vasp_util(
        ["pd", "-d", "Va_O1_0"])
    expected = Namespace(
        dir=Path("Va_O1_0"),
        func=parsed_args.func)
    assert parsed_args == expected


def test_refine_defect_poscar(mocker, simple_cubic):
    mock_structure = mocker.patch(
        "pydefect.cli.vasp.main_vasp_util.Structure")
    mock_loadfn = mocker.patch("pydefect.cli.vasp.main_vasp_util.loadfn")

    parsed_args = parse_args_main_vasp_util(
        ["rdp", "-p", "POSCAR", "-d", "defect_entry.json", "-n", "refined"])
    expected = Namespace(
        structure=mock_structure.from_file.return_value,
        defect_entry=mock_loadfn.return_value,
        poscar_name="refined",
        func=parsed_args.func)
    assert parsed_args == expected
    mock_structure.from_file.assert_called_once_with("POSCAR")
    mock_loadfn.assert_called_once_with("defect_entry.json")


def test_calc_grids(mocker):
    mock_chgcar = mocker.patch("pydefect.cli.vasp.main_vasp_util.Chgcar")
    parsed_args = parse_args_main_vasp_util(
        ["cg", "-c", "CHG"])
    expected = Namespace(
        chgcar=mock_chgcar.from_file.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_chgcar.from_file.assert_called_once_with("CHG")


def test_calc_defect_charge_info(mocker):
    mock_grids = mocker.patch("pydefect.cli.vasp.main_vasp_util.Grids")
    parsed_args = parse_args_main_vasp_util(
        ["cdc", "-p", "PARCHG.0001.ALLK", "-b", "0.3", "-g", "Grids.npz"])
    expected = Namespace(
        parchgs=["PARCHG.0001.ALLK"],
        bin_interval=0.3,
        grids=mock_grids.from_file.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_grids.from_file.assert_called_once_with("Grids.npz")


