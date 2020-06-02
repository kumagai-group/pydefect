# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import numpy as np
from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag
from pydefect.cli.vasp.main_function import make_supercell, make_defect_set, \
    make_defect_entries, make_unitcell, make_competing_phase_dirs, \
    make_chem_pot_diag, make_calc_results, print_file, \
    make_efnv_correction_from_vasp, make_defect_formation_energy
from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection
from pydefect.defaults import defaults
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pymatgen import IStructure, Composition, Structure, Lattice, Element
from pymatgen.io.vasp import Vasprun, Outcar


def test_print():
    args = Namespace(obj="a")
    print_file(args)


def test_make_unitcell(mocker):
    vasprun_band_mock = mocker.Mock(spec=Vasprun, autospec=True)
    outcar_band_mock = mocker.Mock(spec=Outcar, autospec=True)
    outcar_dielectric_mock = mocker.Mock(spec=Outcar, autospec=True)
    args = Namespace(vasprun_band=vasprun_band_mock,
                     outcar_band=outcar_band_mock,
                     outcar_dielectric=outcar_dielectric_mock)
    mock = mocker.patch("pydefect.cli.vasp.main_function.make_unitcell_from_vasp")
    mock.return_value = Unitcell(vbm=1.0,
                                 cbm=2.0,
                                 ele_dielectric_const=np.eye(3),
                                 ion_dielectric_const=np.eye(3))
    make_unitcell(args)
    mock.assert_called_once_with(vasprun_band=vasprun_band_mock,
                                 outcar_band=outcar_band_mock,
                                 outcar_dielectric=outcar_dielectric_mock)


def test_make_competing_phase_dirs(mocker):
    args = Namespace(elements=["Mg", "O"],
                     e_above_hull=0.1)
    mock = mocker.patch("pydefect.cli.vasp.main_function.MpQuery")
    mock_make = mocker.patch("pydefect.cli.vasp.main_function.make_poscars_from_query")
    make_competing_phase_dirs(args)
    mock.assert_called_once_with(element_list=args.elements,
                                 e_above_hull=args.e_above_hull)
    mock_make.assert_called_once_with(materials_query=mock.return_value.materials, path=Path.cwd())


def test_make_chem_pot_diag(mocker, tmpdir):
    def side_effect(key):
        mock_vasprun = mocker.Mock()
        if key == Path("Mg") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("Mg2")
            mock_vasprun.final_energy = -10
        elif key == Path("O") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("O2")
            mock_vasprun.final_energy = -20
        elif key == Path("MgO") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("MgO")
            mock_vasprun.final_energy = -30
        elif key == Path("Al") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("Al")
            mock_vasprun.final_energy = 0
        elif key == Path("MgAl2O4") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("MgAl2O4")
            mock_vasprun.final_energy = -70
        else:
            raise ValueError
        return mock_vasprun

    tmpdir.chdir()
    mock = mocker.patch("pydefect.cli.vasp.main_function.Vasprun", side_effect=side_effect)
    args = Namespace(dirs=[Path("Mg"), Path("MgO"), Path("O")],
                     target=Composition("MgO"))
    make_chem_pot_diag(args)

    args_2 = Namespace(dirs=[Path("Mg"), Path("MgO"), Path("O"), Path("Al"), Path("MgAl2O4")],
                       target=Composition("MgO"))
    make_chem_pot_diag(args_2)


def test_make_supercell_from_matrix(simple_cubic, simple_cubic_2x1x1, tmpdir):
    matrix = [2, 1, 1]
    args = Namespace(unitcell=simple_cubic, matrix=matrix, min_num_atoms=None, max_num_atoms=None)

    tmpdir.chdir()
    make_supercell(args)
    info = loadfn("supercell_info.json")
    assert IStructure.from_file("SPOSCAR") == simple_cubic_2x1x1
    assert info.structure == simple_cubic_2x1x1
    assert info.transform_matrix == [[2, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_make_recommended_supercell(simple_cubic, simple_cubic_2x2x2, tmpdir):
    args = Namespace(unitcell=simple_cubic, matrix=None, min_num_atoms=8, max_num_atoms=8)

    tmpdir.chdir()
    make_supercell(args)
    info = loadfn("supercell_info.json")
    assert IStructure.from_file("SPOSCAR") == simple_cubic_2x2x2
    assert info.structure == simple_cubic_2x2x2
    assert info.transform_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]


def test_make_defect_set(tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    args = Namespace(oxi_states={"He": 1}, dopants=["Li"], kwargs=["Li_H1", "Va_He1", "Va_H1_-1"])
    make_defect_set(args)

    simple_defects = {SimpleDefect(None, "He1", [-1, 0, 1]),
                      SimpleDefect(None, "H1", [-1]),
                      SimpleDefect("Li", "H1", [0]),
                      }
    assert DefectSet.from_yaml() == DefectSet(defects=simple_defects)


def test_make_defect_entries(tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    defect_set = DefectSet({SimpleDefect(None, "He1", [-1, 0])})
    defect_set.to_yaml()
    args = Namespace()
    make_defect_entries(args)
    names = {str(name) for name in Path(".").glob("*")}
    assert names == {'Va_He1_-1', 'defect_in.yaml', 'perfect', 'Va_He1_0', 'supercell_info.json'}

    perfect_structure = Structure.from_file(Path("perfect") / "POSCAR")
    assert perfect_structure == supercell_info.structure

    file_names = {str(file_name.name) for file_name in Path("Va_He1_-1").glob("*")}
    assert file_names == {"POSCAR", "defect_entry.json", "prior_info.yaml"}

    expected = """charge: -1
"""
    assert Path("Va_He1_-1/prior_info.yaml").read_text() == expected


def test_make_calc_results(tmpdir, mocker):
    tmpdir.chdir()
    mock = mocker.patch("pydefect.cli.vasp.main_function.make_calc_results_from_vasp")
    mock_vasprun = mocker.patch("pydefect.cli.vasp.main_function.Vasprun")
    mock_outcar = mocker.patch("pydefect.cli.vasp.main_function.Outcar")
    mock_calc_results = mocker.Mock(spec=CalcResults)
    mock.return_value = mock_calc_results
    args = Namespace(dirs=[Path("a")])
    make_calc_results(args)

    mock_vasprun.assert_called_with(Path("a") / defaults.vasprun)
    mock_outcar.assert_called_with(Path("a") / defaults.outcar)
    mock.assert_called_with(vasprun=mock_vasprun.return_value, outcar=mock_outcar.return_value)
    mock_calc_results.to_json_file.assert_called_with(filename=Path("a") / "calc_results.json")


def test_make_efnv_correction_from_vasp(tmpdir, mocker):
    mock_defect_entry = mocker.Mock(spec=DefectEntry, autospec=True)
    mock_calc_results = mocker.Mock(spec=CalcResults, autospec=True)

    def side_effect(key):
        if str(key) == "Va_O1_2/defect_entry.json":
            mock_defect_entry.full_name = "Va_O1_2"
            mock_defect_entry.charge = 2
            return mock_defect_entry
        elif str(key) == "Va_O1_2/calc_results.json":
            return mock_calc_results
        else:
            raise ValueError

    mock_perfect_calc_results = mocker.Mock(spec=CalcResults)
    mock_loadfn = mocker.patch("pydefect.cli.vasp.main_function.loadfn", side_effect=side_effect)
    mock_unitcell = mocker.Mock(spec=Unitcell)
    mock_make_efnv = mocker.patch("pydefect.cli.vasp.main_function.make_efnv_correction")
    mock_efnv = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    mock_make_efnv.return_value = mock_efnv

    mock_pot_plotter = mocker.patch("pydefect.cli.vasp.main_function.SitePotentialPlotter")

    args = Namespace(dirs=[Path("Va_O1_2")],
                     perfect_calc_results=mock_perfect_calc_results,
                     unitcell=mock_unitcell)

    make_efnv_correction_from_vasp(args)
    mock_loadfn.assert_any_call(Path("Va_O1_2") / "defect_entry.json")
    mock_loadfn.assert_any_call(Path("Va_O1_2") / "calc_results.json")

    mock_make_efnv.assert_called_with(mock_defect_entry.charge, mock_calc_results, mock_perfect_calc_results, mock_unitcell.dielectric_constant)
    mock_efnv.to_json_file.assert_called_with(Path("Va_O1_2") / "correction.json")

    mock_pot_plotter.assert_called_with("Va_O1_2", mock_efnv)
    mock_pot_plotter.return_value.construct_plot.assert_called_once_with()
    mock_pot_plotter.return_value.plt.savefig.assert_called_once_with(fname=Path("Va_O1_2") / "correction.pdf")


def test_make_defect_formation_energy(tmpdir, mocker):
    mock_perfect_calc_results = mocker.Mock(spec=CalcResults)
    mock_perfect_calc_results.structure = Structure(Lattice.cubic(1), species=["H"] * 2, coords=[[0]*3] * 2)
    mock_perfect_calc_results.energy = 10
    mock_perfect_calc_results.vbm = 10
    mock_perfect_calc_results.cbm = 20

    mock_chem_pot_diag = mocker.Mock(spec=ChemPotDiag)
    mock_chem_pot_diag.abs_chem_pot_dict.return_value = {Element.H: 0}

    mock_defect_entry = mocker.Mock(spec=DefectEntry, autospec=True)
    mock_calc_results = mocker.Mock(spec=CalcResults, autospec=True)

    def side_effect(key):
        if str(key) == "Va_O1_2/defect_entry.json":
            mock_defect_entry.name = "Va_H1"
            mock_defect_entry.charge = 2
            return mock_defect_entry
        elif str(key) == "Va_O1_2/calc_results.json":
            mock_calc_results.structure = Structure(Lattice.cubic(1), species=["H"], coords=[[0]*3])
            mock_calc_results.energy = 10
            return mock_calc_results
        elif str(key) == "Va_O1_2/correction.json":
            mock_calc_results.correction_energy = 20
            return mock_calc_results
        else:
            raise ValueError

    mock_loadfn = mocker.patch("pydefect.cli.vasp.main_function.loadfn", side_effect=side_effect)

    mock_unitcell = mocker.Mock(spec=Unitcell)
    mock_unitcell.vbm = 11
    mock_unitcell.cbm = 19

    args = Namespace(dirs=[Path("Va_O1_2")],
                     perfect_calc_results=mock_perfect_calc_results,
                     unitcell=mock_unitcell,
                     chem_pot_diag=mock_chem_pot_diag,
                     label="A",
                     y_range=[-100, 100])
    make_defect_formation_energy(args)

    mock_chem_pot_diag.abs_chem_pot_dict.assert_called_once_with("A")
    mock_loadfn.assert_any_call(Path("Va_O1_2") / "defect_entry.json")
    mock_loadfn.assert_any_call(Path("Va_O1_2") / "calc_results.json")




