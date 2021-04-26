# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.main_util import parse_args_main_util
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection


def test_make_cpd_options():
    parsed_args = parse_args_main_util(["cefm",
                                        "-e", "Mg", "O",
                                        "-a", "atom_energies.yaml"])
    expected = Namespace(
        elements=["Mg", "O"],
        atom_energy_yaml="atom_energies.yaml",
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

