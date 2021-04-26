# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.cli.vasp.utils.main_vasp_util import parse_args_main_vasp_util


def test_make_parchg_dir():
    parsed_args = parse_args_main_vasp_util(
        ["pd", "-d", "Va_O1_0"])
    expected = Namespace(
        dir=Path("Va_O1_0"),
        func=parsed_args.func)
    assert parsed_args == expected


def test_calc_defect_charge_info():
    parsed_args = parse_args_main_vasp_util(
        ["cdc", "-p", "a/PARCHG.0001.ALLK", "-v", "x", "-b", "0.3", "-g", "y"])
    expected = Namespace(
        parchgs=["a/PARCHG.0001.ALLK"],
        vesta_file=Path("x"),
        bin_interval=0.3,
        grids_dirname=Path("y"),
        func=parsed_args.func)
    assert parsed_args == expected


