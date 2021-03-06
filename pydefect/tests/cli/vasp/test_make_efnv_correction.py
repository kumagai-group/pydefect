# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from pydefect.analyzer.calc_results import CalcResults
from pydefect.cli.vasp.make_efnv_correction import \
    calc_max_sphere_radius, make_efnv_correction
from pydefect.corrections.efnv_correction import PotentialSite
from pymatgen.core import IStructure, Lattice


def test_make_efnv_correction(mocker):
    mock_perfect = mocker.Mock(spec=CalcResults, autospec=True)
    mock_defect = mocker.Mock(spec=CalcResults, autospec=True)

    mock_perfect.structure = IStructure(
        Lattice.cubic(10), species=["H"] + ["He"] * 3 + ["Li"],
        coords=[[0, 0, 0], [1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2], [0, 0, 1/2]])
    mock_defect.structure = IStructure(
        Lattice.cubic(10), species=["He"] * 3 + ["Li"],
        coords=[[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2], [0, 0, 1/2]])
    mock_perfect.potentials = [3.0, 4.0, 5.0, 6.0, 7.0]
    mock_defect.potentials = [14.0, 25.0, 36.0, 47.0]

    mock_ewald = mocker.patch("pydefect.cli.vasp.make_efnv_correction.Ewald")
    ewald = mocker.Mock()
    ewald.lattice_energy = 1e3
    ewald.atomic_site_potential.return_value = 1e4
    mock_ewald.return_value = ewald

    efnvc = make_efnv_correction(charge=2,
                                 calc_results=mock_defect,
                                 perfect_calc_results=mock_perfect,
                                 dielectric_tensor=np.eye(3))

    unit_conversion = 180.95128169876497

    assert efnvc.charge == 2
    assert efnvc.point_charge_correction == -4e3 * unit_conversion
    assert efnvc.defect_region_radius == 5.0
    assert efnvc.sites == [PotentialSite("He", 5 * np.sqrt(2), 10.0, 2e4 * unit_conversion),
                           PotentialSite("He", 5 * np.sqrt(2), 20.0, 2e4 * unit_conversion),
                           PotentialSite("He", 5 * np.sqrt(2), 30.0, 2e4 * unit_conversion),
                           PotentialSite("Li", 5.0, 40.0, None),
                           ]


def test_calc_max_sphere_radius():
    lattice_vectors_1 = np.array([[5, 0, 0], [0, 10, 0], [0, 0, 20]])
    lattice_vectors_2 = np.array([[10, 0, 0], [0, 10, 0], [10, 10, 10]])
    radius_1 = calc_max_sphere_radius(lattice_vectors_1)
    radius_2 = calc_max_sphere_radius(lattice_vectors_2)
    assert radius_1 == 10.0
    assert radius_2 == 5.0

