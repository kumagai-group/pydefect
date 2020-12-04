# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag
from pydefect.cli.vasp.main import parse_args
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection
from pydefect.defaults import defaults
from pydefect.input_maker.supercell_info import SupercellInfo
from pymatgen import Composition


def test_print(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["p", "-f", "a.json"])
    expected = Namespace(obj=mock.return_value, func=parsed_args.func)
    assert parsed_args == expected


def test_unitcell(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.Vasprun")
    mock_outcar = mocker.patch("pydefect.cli.vasp.main.Outcar")

    parsed_args = parse_args(["u",
                              "-vb", "vasprun.xml",
                              "-ob", "OUTCAR-1",
                              "-odc", "OUTCAR-2",
                              "-odi", "OUTCAR-3"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        vasprun_band=mock.return_value,
        outcar_band=mock_outcar.return_value,
        outcar_dielectric_clamped=mock_outcar.return_value,
        outcar_dielectric_ionic=mock_outcar.return_value,
        func=parsed_args.func,
    )
    assert parsed_args == expected
    mock.assert_called_once_with("vasprun.xml")
    mock_outcar.assert_any_call("OUTCAR-1")
    mock_outcar.assert_any_call("OUTCAR-2")
    mock_outcar.assert_called_with("OUTCAR-3")


def test_make_poscars_wo_options():
    parsed_args = parse_args(["mp",
                              "-e", "Mg", "O"])
    expected = Namespace(
        elements=["Mg", "O"],
        e_above_hull=defaults.e_above_hull,
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_poscars_w_options():
    parsed_args = parse_args(["mp",
                              "-e", "Mg", "O",
                              "--e_above_hull", "0.1"])
    expected = Namespace(
        elements=["Mg", "O"],
        e_above_hull=0.1,
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_cpd_options():
    parsed_args = parse_args(["mcpd",
                              "-d", "Mg", "O",
                              "-t", "MgO"])
    expected = Namespace(
        dirs=[Path("Mg"), Path("O")],
        target=Composition("MgO"),
        elements=None,
        atom_energy_yaml=None,
        yaml="cpd.yaml",
        update=False,
        func=parsed_args.func)
    assert parsed_args == expected


def test_plot_cpd_options():
    parsed_args = parse_args(["pcpd",
                              "-y", "a.yaml"])
    expected = Namespace(
        yaml="a.yaml",
        func=parsed_args.func)
    assert parsed_args == expected


def test_make_supercell_wo_options(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.IStructure")
    parsed_args = parse_args(["s"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        unitcell=mock.from_file.return_value,
        matrix=None,
        min_num_atoms=50,
        max_num_atoms=300,
        func=parsed_args.func,
    )
    assert parsed_args == expected
    mock.from_file.assert_called_once_with("POSCAR")


def test_make_supercell_w_options(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.IStructure")
    parsed_args = parse_args(["s",
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
        func=parsed_args.func,
    )
    assert parsed_args == expected
    mock.from_file.assert_called_once_with("POSCAR-tmp")


def test_defect_set_wo_options():
    parsed_args = parse_args(["ds"])
    expected = Namespace(
        oxi_states=None,
        dopants=None,
        kwargs=None,
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_defect_set_w_options():
    parsed_args = parse_args(["ds",
                              "-o", "He", "1",
                              "-d", "Li",
                              "-k", "Li_H1", "Va_H1_0"])
    expected = Namespace(
        oxi_states=["He", 1],
        dopants=["Li"],
        kwargs=["Li_H1", "Va_H1_0"],
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_defect_entries():
    parsed_args = parse_args(["de"])
    expected = Namespace(
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_append_interstitial(mocker):

    mock_loadfn = mocker.patch("pydefect.cli.vasp.main.loadfn")
    mock_supercell_info = mocker.Mock(spec=SupercellInfo, autospec=True)
    mock_loadfn.return_value = mock_supercell_info
    mock_structure = mocker.patch("pydefect.cli.vasp.main.Structure")
    parsed_args = parse_args(["ai",
                              "-s", "supercell_info.json",
                              "-p", "POSCAR",
                              "-c", "0.1", "0.2", "0.3"])
    expected = Namespace(
        supercell_info=mock_supercell_info,
        base_structure=mock_structure.from_file.return_value,
        frac_coords=[0.1, 0.2, 0.3],
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_pop_interstitial(mocker):

    mock_loadfn = mocker.patch("pydefect.cli.vasp.main.loadfn")
    mock_supercell_info = mocker.Mock(spec=SupercellInfo, autospec=True)
    mock_loadfn.return_value = mock_supercell_info
    parsed_args = parse_args(["pi",
                              "-s", "supercell_info.json",
                              "-i", "1000"])
    expected = Namespace(
        supercell_info=mock_supercell_info,
        index=1000,
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_calc_results():
    parsed_args = parse_args(["cr", "-d", "Va_O1_0", "Va_O1_1"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_efnv_correction(mocker):
    mock_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_unitcell = mocker.Mock(spec=Unitcell, autospec=True)

    def side_effect(filename):
        if filename == "perfect/calc_results.json":
            return mock_calc_results
        elif filename == "unitcell.json":
            return mock_unitcell
        else:
            raise ValueError

    mocker.patch("pydefect.cli.vasp.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args(["efnv",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json",
                              "-u", "unitcell.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock_calc_results,
        unitcell=mock_unitcell,
        func=parsed_args.func)
    assert parsed_args == expected


def test_gkfo_correction(mocker):
    mock_i_correction = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    mock_i_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_f_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_unitcell = mocker.Mock(spec=Unitcell, autospec=True)

    def side_effect(filename):
        if filename == "a/correction.json":
            return mock_i_correction
        elif filename == "a/calc_results.json":
            return mock_i_calc_results
        elif filename == "a/absorption/calc_results.json":
            return mock_f_calc_results
        elif filename == "unitcell.json":
            return mock_unitcell
        else:
            raise ValueError

    mocker.patch("pydefect.cli.vasp.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args(["gkfo",
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
        unitcell=mock_unitcell,
        func=parsed_args.func)
    assert parsed_args == expected


def test_defect_eigenvalues(mocker):
    mock_calc_results = mocker.Mock(spec=CalcResults, autospec=True)

    def side_effect(filename):
        if filename == "perfect/calc_results.json":
            return mock_calc_results
        else:
            raise ValueError

    mocker.patch("pydefect.cli.vasp.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args(["eig",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock_calc_results,
        func=parsed_args.func)
    assert parsed_args == expected


def test_band_edge_characters(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["mec",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.assert_any_call("perfect/calc_results.json")


def test_band_edge_states(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["es",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-p", "perfect/edge_characters.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_edge_characters=mock.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.assert_any_call("perfect/edge_characters.json")


def test_defect_formation_energy(mocker):
    mock_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_unitcell = mocker.Mock(spec=Unitcell, autospec=True)
    mock_chem_pot_diag = mocker.Mock(spec=ChemPotDiag, autospec=True)

    def side_effect(filename):
        if filename == "perfect/calc_results.json":
            return mock_calc_results
        elif filename == "unitcell.json":
            return mock_unitcell
        else:
            raise ValueError

    mocker.patch("pydefect.cli.vasp.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args(["e",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json",
                              "-u", "unitcell.json",
                              "-c", "cpd.yaml",
                              "-l", "A",
                              "-y", "-5", "5",
                              "-s",
                              "-p"
                              ])

    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock_calc_results,
        unitcell=mock_unitcell,
        cpd_yaml="cpd.yaml",
        label="A",
        y_range=[-5, 5],
        skip_shallow=True,
        print=True,
        supercell_edge=False,
        label_line=True,
        add_charges=True,
        web_gui=False,
        port=8050,
        func=parsed_args.func)
    assert parsed_args == expected
