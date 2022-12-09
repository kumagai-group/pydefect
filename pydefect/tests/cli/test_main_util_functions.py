# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import pytest
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergySummary, DefectEnergies, \
    DefectEnergy
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    CompositionEnergy
from pydefect.cli.main_util_functions import composition_energies_from_mp, \
    make_gkfo_correction_from_vasp, add_interstitials_from_local_extrema, \
    make_defect_vesta_file, show_u_values, show_pinning_levels
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection
from pymatgen.core import Composition


def test_composition_energies_from_mp(mocker, tmpdir):
    tmpdir.chdir()
    args = Namespace(elements=["H"], atom_energy_yaml="a.yaml")
    mock = mocker.patch(
        "pydefect.cli.main_util_functions.make_composition_energies_from_mp")
    mock.return_value = CompositionEnergies(
        {Composition("H"): CompositionEnergy(1.0, "mp-1")})
    composition_energies_from_mp(args)
    mock.assert_called_once_with(elements=["H"], atom_energy_yaml="a.yaml")
    actual = Path("composition_energies.yaml").read_text()
    expected = """H1:
  energy: 1.0
  source: mp-1
"""
    assert actual == expected


#TODO: refactor fixture


def test_show_u_values(defect_energy_summary):
    args = Namespace(defect_energy_summary=defect_energy_summary,
                     with_corrections=True,
                     allow_shallow=True)
    show_u_values(args)


def test_show_pinning_levels(defect_energy_summary):
    args = Namespace(defect_energy_summary=defect_energy_summary,
                     label="A",
                     with_corrections=True,
                     allow_shallow=True)
    show_pinning_levels(args)


def test_add_interstitials_from_local_extrema(mocker):
    mock_supercell_info = mocker.Mock()
    mock_local_extrema = mocker.Mock()
    args = Namespace(supercell_info=mock_supercell_info,
                     local_extrema=mock_local_extrema,
                     indices=[1, 2])
    add_interstitials_from_local_extrema(args)
    mock_local_extrema.append_sites_to_supercell_info.assert_called_once_with(
        mock_supercell_info, [1, 2])


def test_make_defect_vesta_file(mocker):
    mock_defect_structure_info = mocker.Mock()

    def side_effect(key):
        if str(key) == "Va_O1_0/defect_structure_info.json":
            return mock_defect_structure_info
        else:
            print(key)
            raise ValueError

    mock_loadfn = mocker.patch("pydefect.cli.main_util_functions.loadfn",
                               side_effect=side_effect)
    mock_loadfn2 = mocker.patch("pydefect.cli.main_functions.loadfn",
                                side_effect=side_effect)
    mock_make_vesta_file = mocker.patch(
        "pydefect.cli.main_util_functions.MakeDefectVestaFile")
    args = Namespace(
        dirs=[Path("Va_O1_0")],
        check_calc_results=True,
        cutoff=1.0,
        min_displace_w_arrows=2.0,
        arrow_factor=3.0,
        title=None)
    make_defect_vesta_file(args)
    mock_make_vesta_file.assert_called_once_with(mock_defect_structure_info,
                                                 1.0, 2.0, 3.0, None)


def test_make_gkfo_correction_from_vasp(tmpdir, mocker):
    mock_i_correction = mocker.Mock(spec=ExtendedFnvCorrection, autospec=True)
    mock_i_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_f_calc_results = mocker.Mock(spec=CalcResults, autospec=True)
    mock_unitcell = mocker.Mock()

    mock_pot_plotter = mocker.patch(
        "pydefect.cli.main_util_functions.SitePotentialMplPlotter")
    mock_make_gkfo = mocker.patch(
        "pydefect.cli.main_util_functions.make_gkfo_correction")

    args = Namespace(
        initial_efnv_correction=mock_i_correction,
        initial_calc_results=mock_i_calc_results,
        final_calc_results=mock_f_calc_results,
        charge_diff=1,
        unitcell=mock_unitcell)

    make_gkfo_correction_from_vasp(args)
    mock_make_gkfo.assert_called_with(
        efnv_correction=mock_i_correction,
        additional_charge=1,
        final_calc_results=mock_f_calc_results,
        initial_calc_results=mock_i_calc_results,
        diele_tensor=mock_unitcell.dielectric_constant,
        ion_clamped_diele_tensor=mock_unitcell.ele_dielectric_const)


