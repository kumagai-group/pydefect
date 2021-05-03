# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState, \
    BandEdgeOrbitalInfos
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    CompositionEnergy, RelativeEnergies, ChemPotDiag
from pydefect.cli.main_functions import make_standard_and_relative_energies, \
    make_cpd_and_vertices, plot_chem_pot_diag, make_supercell, \
    append_interstitial_to_supercell_info, \
    pop_interstitial_from_supercell_info, make_defect_set, \
    make_band_edge_states_main_func, make_efnv_correction_main_func, \
    calc_defect_structure_info
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pymatgen import IStructure
from pymatgen.core import Composition


def test_make_standard_and_relative_energies(mocker, tmpdir):
    tmpdir.chdir()
    args = Namespace(composition_energies_yaml="a.yaml")
    mock = mocker.patch("pydefect.cli.main_functions.CompositionEnergies")
    mock.from_yaml.return_value = CompositionEnergies(
        {Composition("H"): CompositionEnergy(-1.0),
         Composition("He"): CompositionEnergy(-2.0),
         Composition("HHe"): CompositionEnergy(-5.0)})
    make_standard_and_relative_energies(args)
    mock.from_yaml.assert_called_once_with("a.yaml")
    actual = Path("standard_energies.yaml").read_text()
    expected = """H: -1.0
He: -2.0
"""
    assert actual == expected

    actual = Path("relative_energies.yaml").read_text()
    expected = """HHe: -1.0
"""
    assert actual == expected


def test_make_cpd_and_vertices(mocker, tmpdir):
    tmpdir.chdir()
    mock = mocker.patch("pydefect.cli.main_functions.RelativeEnergies")
    mock.from_yaml.return_value = RelativeEnergies({"HHe": -1.0})

    args = Namespace(rel_energy_yaml="a.yaml", target="HHe", elements=None)
    make_cpd_and_vertices(args)

    actual = Path("target_vertices.yaml").read_text()
    expected = """target: HHe
A:
  chem_pot:
    H: 0.0
    He: -2.0
  competing_phases:
  - H
  impurity_phases: []
B:
  chem_pot:
    H: -2.0
    He: 0.0
  competing_phases:
  - He
  impurity_phases: []
"""
    assert actual == expected


def test_plot_chem_pot_diag(tmpdir):
    tmpdir.chdir()
    chem_pot_diag = ChemPotDiag(vertex_elements=["Mg", "O"],
                                polygons={"Mg": [[0.0, -3.0], [0.0, -1.0]],
                                          "O": [[-3.0, 0.0], [-1.0, 0.0]],
                                          "MgO": [[-1.0, 0.0], [0.0, -1.0]]})
    args = Namespace(chem_pot_diag=chem_pot_diag)
    plot_chem_pot_diag(args)


def test_make_supercell_from_matrix(simple_cubic, simple_cubic_2x1x1, tmpdir):
    tmpdir.chdir()
    matrix = [2, 1, 1]
    args = Namespace(unitcell=simple_cubic, matrix=matrix,
                     min_num_atoms=None, max_num_atoms=None)
    make_supercell(args)
    info = loadfn("supercell_info.json")
    assert IStructure.from_file("SPOSCAR") == simple_cubic_2x1x1
    assert info.structure == simple_cubic_2x1x1
    assert info.transformation_matrix == [[2, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_make_recommended_supercell(simple_cubic, simple_cubic_2x2x2, tmpdir):
    tmpdir.chdir()
    args = Namespace(unitcell=simple_cubic, matrix=None,
                     min_num_atoms=8, max_num_atoms=8)
    make_supercell(args)
    info = loadfn("supercell_info.json")
    assert IStructure.from_file("SPOSCAR") == simple_cubic_2x2x2
    assert info.structure == simple_cubic_2x2x2
    assert info.transformation_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]


def test_add_interstitials(mocker):
    mock_1 = mocker.Mock()
    mock_2 = mocker.Mock()
    mock_3 = mocker.Mock()
    mock_4 = mocker.Mock()
    args = Namespace(supercell_info=mock_1, base_structure=mock_2,
                     frac_coords=mock_3, info=mock_4)
    mock = mocker.patch("pydefect.cli.main_functions.append_interstitial")

    append_interstitial_to_supercell_info(args)
    mock.assert_called_once_with(mock_1, mock_2, [mock_3], mock_4)
    mock.return_value.to_json_file.assert_called_once_with()


def test_pop_interstitials(mocker):
    mock_si = mocker.MagicMock()
    args = Namespace(supercell_info=mock_si, index=1000, pop_all=False)

    pop_interstitial_from_supercell_info(args)
    assert isinstance(mock_si.interstitials, object)
    mock_si.interstitials.pop.assert_called_once_with(999)
    mock_si.to_json_file.assert_called_once_with()


def test_pop_interstitials_pop_all(mocker):
    mock_si = mocker.MagicMock()
    mock_si.interstitials = ["a", "b"]
    args = Namespace(supercell_info=mock_si, index=None, pop_all=True)

    pop_interstitial_from_supercell_info(args)
    assert isinstance(mock_si.interstitials, object)
    mock_si.interstitials = []


@pytest.mark.parametrize("oxi_states,he_vacancy_charge",
                         ([None, [0]], [["He", 1, "Li", 1], [-1, 0, 1]]))
def test_make_defect_set(oxi_states, he_vacancy_charge, tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    args = Namespace(oxi_states=oxi_states, dopants=["Li"],
                     kwargs=["Li_H1", "Va_He1", "Va_H1_-1"])
    make_defect_set(args)

    simple_defects = {SimpleDefect(None, "He1", he_vacancy_charge),
                      SimpleDefect(None, "H1", [-1]),
                      SimpleDefect("Li", "H1", [0])}

    DefectSet(defects=simple_defects).to_yaml("expected.yaml")
    actual = Path("defect_in.yaml").read_text()
    expected = Path("expected.yaml").read_text()
    assert actual == expected


def test_calc_defect_structure_info(mocker, tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    mock_calc_results = mocker.Mock()
    mock_defect_entry = mocker.Mock()

    def side_effect(key):
        if str(key) == "Va_O1_2/calc_results.json":
            return mock_calc_results
        elif str(key) == "Va_O1_2/defect_entry.json":
            return mock_defect_entry
        else:
            print(key)
            raise ValueError

    mock_loadfn = mocker.patch(
        "pydefect.cli.main_functions.loadfn", side_effect=side_effect)
    mock_structure_info = mocker.patch(
        "pydefect.cli.main_functions.MakeDefectStructureInfo")
    info = loadfn("supercell_info.json")

    args = Namespace(supercell_info=info, check_calc_results=True,
                     dirs=[Path("Va_O1_2")], dist_tolerance=0.1, symprec=0.2)
    calc_defect_structure_info(args)
    mock_structure_info.assert_called_with(
        supercell_info.structure,
        mock_defect_entry.structure,
        mock_calc_results.structure,
        dist_tol=0.1,
        symprec=0.2)
    mock_structure_info.return_value.defect_structure_info.to_json_file.assert_called_once_with(
        "Va_O1_2/defect_structure_info.json")


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
    mock_loadfn = mocker.patch(
        "pydefect.cli.main_functions.loadfn", side_effect=side_effect)
    mock_unitcell = mocker.Mock(spec=Unitcell)
    mock_make_efnv = mocker.patch(
        "pydefect.cli.main_functions.make_efnv_correction")
    mock_efnv = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    mock_make_efnv.return_value = mock_efnv

    mock_site_pot_plotter = mocker.patch(
        "pydefect.cli.main_functions.SitePotentialMplPlotter")
    plotter = mock_site_pot_plotter.from_efnv_corr.return_value
    args = Namespace(dirs=[Path("Va_O1_2")],
                     check_calc_results=True,
                     perfect_calc_results=mock_perfect_calc_results,
                     unitcell=mock_unitcell)

    make_efnv_correction_main_func(args)
    mock_loadfn.assert_any_call(Path("Va_O1_2") / "defect_entry.json")
    mock_loadfn.assert_any_call(Path("Va_O1_2") / "calc_results.json")

    mock_make_efnv.assert_called_with(
        mock_defect_entry.charge, mock_calc_results, mock_perfect_calc_results,
        mock_unitcell.dielectric_constant)
    mock_efnv.to_json_file.assert_called_with(
        Path("Va_O1_2") / "correction.json")

    mock_site_pot_plotter.from_efnv_corr.assert_called_with(
        title="Va_O1_2", efnv_correction=mock_efnv)
    plotter.construct_plot.assert_called_once_with()


def test_make_band_edge_states(mocker):
    mock_perfect_edge_states = mocker.Mock(
        spec=PerfectBandEdgeState, autospec=True)
    mocker_band_edge_orbital_infos = mocker.Mock(spec=BandEdgeOrbitalInfos,
                                                 autospec=True)

    def side_effect(key):
        if str(key) == "Va_O1_2/band_edge_orbital_infos.json":
            return mocker_band_edge_orbital_infos
        else:
            raise FileNotFoundError

    mocker_loadfn = mocker.patch(
        "pydefect.cli.main_functions.loadfn", side_effect=side_effect)
    mocker_make_edge_states = mocker.patch(
        "pydefect.cli.main_functions.make_band_edge_states")

    args = Namespace(dirs=[Path("Va_O1_2")],
                     check_calc_results=True,
                     p_state=mock_perfect_edge_states)
    make_band_edge_states_main_func(args)
    mocker_make_edge_states.assert_called_with(mocker_band_edge_orbital_infos,
                                               mock_perfect_edge_states,
                                               None)

