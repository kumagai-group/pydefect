# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.cli.vasp.main import parse_args
from pydefect.defaults import defaults
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
                              "-od", "OUTCAR-2"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        vasprun_band=mock.return_value,
        outcar_band=mock_outcar.return_value,
        outcar_dielectric=mock_outcar.return_value,
        func=parsed_args.func,
    )
    assert parsed_args == expected
    mock.assert_called_once_with("vasprun.xml")
    mock_outcar.assert_any_call("OUTCAR-1")
    mock_outcar.assert_called_with("OUTCAR-2")


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


def test_cpd_wo_options():
    parsed_args = parse_args(["cpd",
                              "-d", "Mg", "O"])
    expected = Namespace(
        dirs=[Path("Mg"), Path("O")],
        target=None,
        func=parsed_args.func)
    assert parsed_args == expected


def test_cpd_w_options():
    parsed_args = parse_args(["cpd",
                              "-d", "Mg", "O",
                              "-t", "MgO"])
    expected = Namespace(
        dirs=[Path("Mg"), Path("O")],
        target=Composition("MgO"),
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
                              "--min_num_atoms", "1000",
                              "--max_num_atoms", "2000"])
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
        oxi_states=["He", "1"],
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


def test_calc_results():
    parsed_args = parse_args(["cr", "-d", "Va_O1_0", "Va_O1_1"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        func=parsed_args.func,
    )
    assert parsed_args == expected


def test_efnv_correction(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["efnv",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json",
                              "-u", "unitcell.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock.return_value,
        unitcell=mock.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.assert_any_call("perfect/calc_results.json")
    mock.assert_any_call("unitcell.json")


def test_defect_eigenvalues(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["eig",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json",
                              "-u", "unitcell.json"])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock.return_value,
        unitcell=mock.return_value,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.assert_any_call("perfect/calc_results.json")
    mock.assert_any_call("unitcell.json")


def test_band_edge_characters(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["ec",
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
    mock = mocker.patch("pydefect.cli.vasp.main.loadfn")
    parsed_args = parse_args(["e",
                              "-d", "Va_O1_0", "Va_O1_1",
                              "-pcr", "perfect/calc_results.json",
                              "-u", "unitcell.json",
                              "-c", "chem_pot_diag.json",
                              "-l", "A",
                              "-y", "-5", "5",
                              "-s",
                              ])
    expected = Namespace(
        dirs=[Path("Va_O1_0"), Path("Va_O1_1")],
        perfect_calc_results=mock.return_value,
        unitcell=mock.return_value,
        chem_pot_diag=mock.return_value,
        label="A",
        y_range=[-5, 5],
        skip_shallow=True,
        func=parsed_args.func)
    assert parsed_args == expected
    mock.assert_any_call("perfect/calc_results.json")
    mock.assert_any_call("unitcell.json")
    mock.assert_any_call("chem_pot_diag.json")
