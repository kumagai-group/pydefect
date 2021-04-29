# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.calc_results import CalcResults
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    CompositionEnergy
from pydefect.cli.main_util_functions import composition_energies_from_mp, \
    make_gkfo_correction_from_vasp, add_interstitials_from_local_extrema
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


def test_add_interstitials_from_local_extrema(mocker):
    mock_supercell_info = mocker.Mock()
    mock_local_extrema = mocker.Mock()
    args = Namespace(supercell_info=mock_supercell_info,
                     local_extrema=mock_local_extrema,
                     indices=[1, 2])
    add_interstitials_from_local_extrema(args)
    mock_local_extrema.append_sites_to_supercell_info.assert_called_once_with(
        mock_supercell_info, [1, 2])


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


