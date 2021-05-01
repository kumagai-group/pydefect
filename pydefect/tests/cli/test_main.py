# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergySummary
from pydefect.chem_pot_diag.chem_pot_diag import TargetVertices
from pydefect.cli.main import parse_args_main
from pydefect.input_maker.supercell_info import SupercellInfo


def test_make_standard_and_relative_energies_options():
    parsed_args = parse_args_main(["sre"])
    expected = Namespace(
        composition_energies_yaml="composition_energies.yaml",
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_cpd_and_vertices_options():
    parsed_args = parse_args_main(["cv",
                                   "-t", "MgO",
                                   "-e", "Mg", "O", "Al"])
    expected = Namespace(
        rel_energy_yaml="relative_energies.yaml",
        target="MgO",
        elements=["Mg", "O", "Al"],
        func=parsed_args.func)
    assert parsed_args == expected


def test_plot_cpd_options(mocker):
    mock_loadfn = mocker.patch("pydefect.cli.main.loadfn")
    parsed_args = parse_args_main(["pc"])
    expected = Namespace(
        chem_pot_diag=mock_loadfn.return_value,
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_supercell_wo_options(mocker):
    mock = mocker.patch("pydefect.cli.main.IStructure")
    parsed_args = parse_args_main(["s", "-p", "POSCAR-tmp"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        unitcell=mock.from_file.return_value,
        matrix=None,
        min_num_atoms=50,
        max_num_atoms=300,
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_supercell_w_options(mocker):
    mock = mocker.patch("pydefect.cli.main.IStructure")
    parsed_args = parse_args_main(["s",
                                   "-p", "POSCAR-tmp",
                                   "--matrix", "1", "2", "3",
                                   "--min_atoms", "1000",
                                   "--max_atoms", "2000"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        unitcell=mock.from_file.return_value,
        matrix=[1, 2, 3],
        min_num_atoms=1000,
        max_num_atoms=2000,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.from_file.assert_called_once_with("POSCAR-tmp")


def test_append_interstitial(mocker):

    mock_loadfn = mocker.patch("pydefect.cli.main.loadfn")
    mock_supercell_info = mocker.Mock(spec=SupercellInfo, autospec=True)
    mock_loadfn.return_value = mock_supercell_info
    mock_structure = mocker.patch("pydefect.cli.main.Structure")
    parsed_args = parse_args_main(["ai",
                                   "-s", "supercell_info.json",
                                   "-p", "POSCAR",
                                   "-c", "0.1", "0.2", "0.3",
                                   "-i", "info"])
    expected = Namespace(
        supercell_info=mock_supercell_info,
        base_structure=mock_structure.from_file.return_value,
        frac_coords=[0.1, 0.2, 0.3],
        info="info",
        func=parsed_args.func)
    assert parsed_args == expected


def test_pop_interstitial(mocker):

    mock_loadfn = mocker.patch("pydefect.cli.main.loadfn")
    mock_supercell_info = mocker.Mock(spec=SupercellInfo, autospec=True)
    mock_loadfn.return_value = mock_supercell_info
    parsed_args = parse_args_main(["pi",
                                   "-s", "supercell_info.json",
                                   "-i", "1000",
                                   "--pop_all"])
    expected = Namespace(
        supercell_info=mock_supercell_info,
        index=1000,
        pop_all=True,
        func=parsed_args.func)
    assert parsed_args == expected


def test_defect_set_wo_options():
    parsed_args = parse_args_main(["ds"])
    expected = Namespace(
        oxi_states=None,
        dopants=None,
        kwargs=None,
        func=parsed_args.func)
    assert parsed_args == expected


def test_defect_set_w_options():
    parsed_args = parse_args_main(["ds",
                                   "-o", "He", "1",
                                   "-d", "Li",
                                   "-k", "Li_H1", "Va_H1_0"])
    expected = Namespace(
        oxi_states=["He", 1],
        dopants=["Li"],
        kwargs=["Li_H1", "Va_H1_0"],
        func=parsed_args.func)
    assert parsed_args == expected


def test_defect_structure_info(mocker):
    mock_supercell_info = mocker.Mock(spec=SupercellInfo, autospec=True)

    def side_effect(filename):
        if filename == "supercell_info.json":
            return mock_supercell_info
        else:
            raise ValueError

    mock_loadfn = mocker.patch("pydefect.cli.main.loadfn",
                               side_effect=side_effect)
    parsed_args = parse_args_main(["dsi",
                                   "-s", "supercell_info.json",
                                   "-d", "Va_O1_0",
                                   "-nccr",
                                   "-dt", "1.0",
                                   "--symprec", "2.0"])
    expected = Namespace(
        supercell_info=mock_supercell_info,
        check_calc_results=False,
        dirs=[Path("Va_O1_0")],
        dist_tolerance=1.0,
        symprec=2.0,
        func=parsed_args.func)
    assert parsed_args == expected


def test_efnv_correction(mocker):
    mock_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")

    def side_effect(filename):
        if filename == "perfect/calc_results.json":
            return mock_calc_results
        else:
            raise ValueError

    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args_main(["efnv",
                                   "-d", "Va_O1_0", "Va_O1_1",
                                   "-nccr",
                                   "-pcr", "perfect/calc_results.json",
                                   "-u", "unitcell.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        check_calc_results=False,
        perfect_calc_results=mock_calc_results,
        unitcell=mock_unitcell.from_yaml.return_value,
        func=parsed_args.func)
    assert parsed_args == expected


def test_band_edge_states(mocker):
    mock = mocker.patch("pydefect.cli.main.loadfn")
    parsed_args = parse_args_main([
        "bes",
        "-d", "Va_O1_0", "Va_O1_1",
        "-nccr",
        "-p", "perfect/perfect_band_edge_state.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        check_calc_results=False,
        p_state=mock.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.assert_any_call("perfect/perfect_band_edge_state.json")


def test_defect_energy_infos(mocker):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_loadfn = mocker.patch("pydefect.cli.main.loadfn")
    mock_std_energy = mocker.patch("pydefect.cli.main.StandardEnergies")
    parsed_args = parse_args_main([
        "dei",
        "-d", "Va_O1_0", "Va_O1_1",
        "-nccr",
        "-u", "unitcell.json",
        "-p", "perfect/perfect_band_edge_state.json",
        "-s", "standard_energies.yaml"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        check_calc_results=False,
        unitcell=mock_unitcell.from_yaml.return_value,
        perfect_calc_results=mock_loadfn.return_value,
        std_energies=mock_std_energy.from_yaml.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_loadfn.assert_called_once_with("perfect/perfect_band_edge_state.json")
    mock_std_energy.from_yaml.assert_called_once_with("standard_energies.yaml")


def test_defect_energy_summary(mocker):
    mock_pbes = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")

    def side_effect(filename):
        if filename == "perfect/calc_results.json":
            return mock_pbes
        else:
            print(filename)
            raise ValueError

    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args_main(["des",
                                   "-d", "Va_O1_0", "Va_O1_1",
                                   "-p", "perfect/calc_results.json",
                                   "-u", "unitcell.json",
                                   "-t", "target_vertices.yaml"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        p_state=mock_pbes,
        unitcell=mock_unitcell.from_yaml.return_value,
        target_vertices_yaml="target_vertices.yaml",
        func=parsed_args.func)
    assert parsed_args == expected


def test_plot_defect_formation_energy(mocker):
    mock_summary = mocker.Mock(spec=DefectEnergySummary, autospec=True)

    def side_effect(filename):
        if filename == "defect_energy_summary.json":
            return mock_summary
        else:
            raise ValueError

    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args_main([
        "pe", "-d", "defect_energy_summary.json", "-l", "A", "-y", "-5", "5",
        "--allow_shallow", "--no_corrections", "--no_label_line",
        "--no_add_charges"])

    expected = Namespace(
        defect_energy_summary=mock_summary,
        label="A",
        y_range=[-5, 5],
        allow_shallow=True,
        with_corrections=False,
        label_line=False,
        add_charges=False,
        func=parsed_args.func)
    assert parsed_args == expected

