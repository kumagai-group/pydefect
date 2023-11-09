# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergySummary
from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.main_util import parse_args_main_util
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection
from pydefect.defaults import defaults


def test_make_cpd_options():
    parsed_args = parse_args_main_util(["cefm",
                                        "-e", "Mg", "O",
                                        "-a", "atom_energies.yaml"])
    expected = Namespace(
        elements=["Mg", "O"],
        atom_energy_yaml="atom_energies.yaml",
        func=parsed_args.func)
    assert parsed_args == expected


def test_show_u_values(mocker):
    mock_summary = mocker.Mock(spec=DefectEnergySummary, autospec=True)

    def side_effect(filename):
        if filename == "defect_energy_summary.json":
            return mock_summary
        else:
            raise ValueError

    # because add_sub_parser is imported from main, patch to loadfn should be
    # done in main.
    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args_main_util(["u",
                                        "-d", "defect_energy_summary.json",
                                        "--allow_shallow",
                                        "--no_corrections",
                                        "--label", "A"])

    expected = Namespace(
        defect_energy_summary=mock_summary,
        allow_shallow=True,
        with_corrections=False,
        label="A",
        func=parsed_args.func)
    assert parsed_args == expected


def test_add_interstitials_from_local_extrema(mocker):

    mock_loadfn_main = mocker.patch("pydefect.cli.main.loadfn")
    mock_loadfn_main_util = mocker.patch("pydefect.cli.main_util.loadfn")

    parsed_args = parse_args_main_util(["ai",
                                        "-s", "supercell_info.json",
                                        "--local_extrema", "local_extrema.json",
                                        "--indices", "1", "2"])
    expected = Namespace(
        supercell_info=mock_loadfn_main.return_value,
        local_extrema=mock_loadfn_main_util.return_value,
        indices=[1, 2],
        func=parsed_args.func)
    assert parsed_args == expected
    mock_loadfn_main.assert_called_once_with("supercell_info.json")
    mock_loadfn_main_util.assert_called_once_with("local_extrema.json")


def test_defect_vesta_file_wo_options():
    parsed_args = parse_args_main_util(["dvf", "-d", "Va_O1_0", "Va_O1_1"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        cutoff=defaults.show_structure_cutoff,
        min_displace_w_arrows=0.1,
        arrow_factor=3.0,
        title=None,
        verbose=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_defect_vesta_file_w_options():
    parsed_args = parse_args_main_util(["dvf",
                                        "-d", "Va_O1_0", "Va_O1_1",
                                        "--cutoff", "1.0",
                                        "--min_displace_w_arrows", "2.0",
                                        "--arrow_factor", "10.0",
                                        "--title", "title"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        cutoff=1.0,
        min_displace_w_arrows=2.0,
        arrow_factor=10.0,
        title="title",
        verbose=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_gkfo_correction(mocker):
    mock_i_correction = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    mock_i_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_f_calc_results = mocker.Mock(spec=CalcResults, autospec=True)

    def side_effect(filename):
        if filename == "a/correction.json":
            return mock_i_correction
        elif filename == "a/calc_results.json":
            return mock_i_calc_results
        elif filename == "a/absorption/calc_results.json":
            return mock_f_calc_results
        else:
            raise ValueError

    mocker.patch("pydefect.cli.main_util.loadfn", side_effect=side_effect)
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")

    parsed_args = parse_args_main_util([
        "gkfo",
        "-iefnv", "a/correction.json",
        "-cd", "1",
        "-icr", "a/calc_results.json",
        "-fcr", "a/absorption/calc_results.json",
        "-u", "unitcell.json"])

    expected = Namespace(
        initial_efnv_correction=mock_i_correction,
        initial_calc_results=mock_i_calc_results,
        final_calc_results=mock_f_calc_results,
        charge_diff=1,
        unitcell=mock_unitcell.from_yaml.return_value,
        func=parsed_args.func)

    assert parsed_args == expected

