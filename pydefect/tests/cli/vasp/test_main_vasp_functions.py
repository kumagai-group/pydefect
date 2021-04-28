# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import numpy as np
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.vasp.main_vasp_functions import make_unitcell, \
    make_competing_phase_dirs, make_composition_energies, make_defect_entries, \
    make_calc_results, make_band_edge_orb_infos_and_eigval_plot, \
    make_perfect_band_edge_state, make_local_extrema
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pymatgen.core import Composition, Structure
from pymatgen.io.vasp import Vasprun, Outcar, VolumetricData
from vise.defaults import defaults


def test_make_unitcell(tmpdir, mocker):
    tmpdir.chdir()
    vasprun_band_mock = mocker.Mock(spec=Vasprun, autospec=True)
    outcar_band_mock = mocker.Mock(spec=Outcar, autospec=True)
    outcar_dielectric_mock = mocker.Mock(spec=Outcar, autospec=True)
    args = Namespace(vasprun_band=vasprun_band_mock,
                     outcar_band=outcar_band_mock,
                     outcar_dielectric_clamped=outcar_dielectric_mock,
                     outcar_dielectric_ionic=outcar_dielectric_mock,
                     name="H")

    mock = mocker.patch(
        "pydefect.cli.vasp.main_vasp_functions.make_unitcell_from_vasp")
    mock.return_value = Unitcell(system="H",
                                 vbm=1.0,
                                 cbm=2.0,
                                 ele_dielectric_const=np.eye(3).tolist(),
                                 ion_dielectric_const=np.eye(3).tolist())
    make_unitcell(args)
    mock.assert_called_once_with(
        vasprun_band=vasprun_band_mock,
        outcar_band=outcar_band_mock,
        outcar_dielectric_clamped=outcar_dielectric_mock,
        outcar_dielectric_ionic=outcar_dielectric_mock,
        system_name="H")


def test_make_competing_phase_dirs(mocker):
    args = Namespace(elements=["Mg", "O"],
                     e_above_hull=0.1)
    mock = mocker.patch("pydefect.cli.vasp.main_vasp_functions.MpQuery")
    mock_make = mocker.patch(
        "pydefect.cli.vasp.main_vasp_functions.make_poscars_from_query")
    make_competing_phase_dirs(args)
    mock.assert_called_once_with(element_list=args.elements,
                                 e_above_hull=args.e_above_hull)
    mock_make.assert_called_once_with(
        materials_query=mock.return_value.materials, path=Path.cwd())


def test_make_composition_energies(mocker, tmpdir):
    def side_effect(key):
        mock_vasprun = mocker.Mock()
        if key == Path("Mg") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("Mg2")
            mock_vasprun.final_energy = -10
        elif key == Path("O") / defaults.vasprun:
            mock_vasprun.final_structure.composition = Composition("O2")
            mock_vasprun.final_energy = -20
        else:
            raise ValueError
        return mock_vasprun

    tmpdir.chdir()
    print(tmpdir)
    mock = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Vasprun",
                        side_effect=side_effect)
    args = Namespace(yaml_file=None, dirs=[Path("Mg"), Path("O")])
    make_composition_energies(args)
    actual = Path("composition_energies.yaml").read_text()
    expected = """Mg2:
  energy: -10.0
  source: local
O2:
  energy: -20.0
  source: local
"""
    assert actual == expected


def test_make_local_extrema(tmpdir, mocker, simple_cubic):
    print(tmpdir)
    tmpdir.chdir()

    mock_params = mocker.patch("pydefect.cli.vasp.main_vasp_functions.VolumetricDataAnalyzeParams")
    mock_make_extrema = mocker.patch("pydefect.cli.vasp.main_vasp_functions.make_local_extrema_from_volumetric_data")
    volumetric_data = VolumetricData(simple_cubic,
                                     data={"total": np.array([[[0.0]]])})
    args = Namespace(volumetric_data=[volumetric_data],
                     find_min=True,
                     info="a",
                     threshold_frac=None,
                     threshold_abs=None,
                     min_dist=0.1,
                     tol=0.2,
                     radius=0.3)
    make_local_extrema(args)
    mock_params.assert_called_once_with(None, None, 0.1, 0.2, 0.3)
    mock_make_extrema.assert_called_once_with(volumetric_data=volumetric_data,
                                              params=mock_params.return_value,
                                              info="a",
                                              find_min=True)
    mock_make_extrema.return_value.to_json_file.assert_called_once_with()


def test_make_defect_entries(tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    defect_set = DefectSet({SimpleDefect(None, "He1", [-1, 0])})
    defect_set.to_yaml()
    args = Namespace()
    make_defect_entries(args)
    names = {str(name) for name in Path(".").glob("*")}
    assert names == {'Va_He1_-1', 'defect_in.yaml', 'perfect', 'Va_He1_0',
                     'supercell_info.json'}

    perfect_structure = Structure.from_file(str(Path("perfect") / "POSCAR"))
    assert perfect_structure == supercell_info.structure

    file_names = {str(file_name.name)
                  for file_name in Path("Va_He1_-1").glob("*")}
    assert file_names == {"POSCAR", "defect_entry.json", "prior_info.yaml"}

    expected = """charge: -1
"""
    assert Path("Va_He1_-1/prior_info.yaml").read_text() == expected


def test_make_calc_results(tmpdir, mocker):
    tmpdir.chdir()
    mock = mocker.patch(
        "pydefect.cli.vasp.main_vasp_functions.make_calc_results_from_vasp")
    mock_vasprun = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Vasprun")
    mock_outcar = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Outcar")
    mock_calc_results = mocker.Mock(spec=CalcResults)
    mock.return_value = mock_calc_results
    args = Namespace(dirs=[Path("a")])
    make_calc_results(args)

    mock_vasprun.assert_called_with(Path("a") / defaults.vasprun)
    mock_outcar.assert_called_with(Path("a") / defaults.outcar)
    mock.assert_called_with(vasprun=mock_vasprun.return_value,
                            outcar=mock_outcar.return_value)
    mock_calc_results.to_json_file.assert_called_with("a/calc_results.json")


def test_make_band_edge_orb_infos_and_eigval_plot(mocker):
    mock_procar = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Procar")
    mock_vasprun = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Vasprun")

    mock_p_state = mocker.Mock()
    mock_p_state.vbm_info.energy = 10
    mock_p_state.cbm_info.energy = 20

    mock_defect_entry = mocker.Mock(spec=DefectEntry, autospec=True)

    mock_structure_info = mocker.Mock(spec=DefectStructureInfo, autospec=True)

    mock_make_orbital_infos = mocker.patch(
        "pydefect.cli.vasp.main_vasp_functions.make_band_edge_orbital_infos")
    mock_eigval_plotter = mocker.patch(
        "pydefect.cli.vasp.main_vasp_functions.EigenvalueMplPlotter")

    def side_effect(key):
        if str(key) == "Va_O1_2/defect_entry.json":
            mock_defect_entry.name = "Va_O1"
            mock_defect_entry.charge = 2
            return mock_defect_entry
        elif str(key) == "Va_O1_2/defect_structure_info.json":
            return mock_structure_info
        else:
            raise ValueError

    mock_loadfn = mocker.patch("pydefect.cli.vasp.main_vasp_functions.loadfn",
                               side_effect=side_effect)

    args = Namespace(dirs=[Path("Va_O1_2")],
                     p_state=mock_p_state)
    make_band_edge_orb_infos_and_eigval_plot(args)

    mock_procar.assert_called_with(Path("Va_O1_2") / defaults.procar)
    mock_vasprun.assert_called_with(Path("Va_O1_2") / defaults.vasprun)

    mock_make_orbital_infos.assert_called_with(
        mock_procar.return_value, mock_vasprun.return_value, 10, 20,
        mock_structure_info)

    mock_loadfn.assert_any_call(Path("Va_O1_2") / "defect_entry.json")

    mock_eigval_plotter.assert_called_with(
        title="Va_O1",
        band_edge_orb_infos=mock_make_orbital_infos.return_value,
        supercell_vbm=10,
        supercell_cbm=20)


def test_make_perfect_band_edge_state(mocker):
    mock_vasprun = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Vasprun")
    mock_procar = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Procar")
    mock_outcar = mocker.patch("pydefect.cli.vasp.main_vasp_functions.Outcar")
    mock_make_perf_be_state = mocker.patch(
        "pydefect.cli.vasp.main_vasp_functions.make_perfect_band_edge_state_from_vasp")

    args = Namespace(dir=Path("perfect"))
    make_perfect_band_edge_state(args)

    mock_vasprun.assert_called_with(Path("perfect") / defaults.vasprun)
    mock_procar.assert_called_with(Path("perfect") / defaults.procar)
    mock_outcar.assert_called_with(Path("perfect") / defaults.outcar)
    mock_make_perf_be_state.assert_called_with(mock_procar.return_value,
                                               mock_vasprun.return_value,
                                               mock_outcar.return_value)
