# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import EdgeInfo, BandEdgeOrbitalInfos, \
    PerfectBandEdgeState, OrbitalInfo, BandEdgeStates, BandEdgeState, \
    LocalizedOrbital
from pydefect.analyzer.defect_charge_info import DefectChargeInfo
from pydefect.analyzer.make_band_edge_states import make_band_edge_states, \
    orbital_diff


@pytest.fixture
def p_edge_state():
    vbm_info = EdgeInfo(band_idx=10, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=OrbitalInfo(energy=-1.0, orbitals={"Mn": [0.5, 0.6, 0.0, 0.0], "O": [0.0, 0.0, 0.0, 0.0]}, occupation=1.0, participation_ratio=0.1))
    cbm_info = EdgeInfo(band_idx=11, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=OrbitalInfo(energy=1.0, orbitals={"Mn": [0.0, 0.0, 0.0, 0.0], "O": [0.1, 0.2, 0.0, 0.0]}, occupation=0.0, participation_ratio=0.1))
    return PerfectBandEdgeState(vbm_info=vbm_info, cbm_info=cbm_info)


@pytest.fixture
def orb_infos():
    orbital_infos = [[[
        OrbitalInfo(energy=-1.1, orbitals={"Mn": [0.5, 0.6, 0.0, 0.0],
                                           "O": [0.0, 0.0, 0.0, 0.0]},
                    occupation=1.0, participation_ratio=0.1),
        OrbitalInfo(energy=-0.9, orbitals={"Mn": [0.5, 0.7, 0.0, 0.0],
                                           "O": [0.0, 0.0, 0.0, 0.0]},
                    occupation=1.0, participation_ratio=0.1),  # vbm
        OrbitalInfo(energy=0.0,  orbitals={"Mn": [0.1, 0.2, 0.0, 0.0],
                                           "O": [0.3, 0.4, 0.0, 0.0]},
                    occupation=1.0, participation_ratio=0.1),  # in-gap
        OrbitalInfo(energy=1.0,  orbitals={"Mn": [0.5, 0.8, 0.0, 0.0],
                                           "O": [0.0, 0.0, 0.0, 0.0]},
                    occupation=0.0, participation_ratio=0.1),  # in-gap
        OrbitalInfo(energy=1.2,  orbitals={"Mn": [0.0, 0.0, 0.0, 0.0],
                                           "O": [0.1, 0.3, 0.0, 0.0]},
                    occupation=0.0, participation_ratio=0.1)]]]  # cbm
    return BandEdgeOrbitalInfos(kpt_coords=[(0.0, 0.0, 0.0)],
                                kpt_weights=[1.0],
                                orbital_infos=orbital_infos,
                                lowest_band_index=8,
                                fermi_level=0.5)


def test_make_band_edge_state(p_edge_state, orb_infos):
    actual = make_band_edge_states(orb_infos, p_edge_state)

    vbm_info = EdgeInfo(band_idx=9, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=OrbitalInfo(
                            energy=-0.9, orbitals={"Mn": [0.5, 0.7, 0.0, 0.0],
                                                   "O": [0.0, 0.0, 0.0, 0.0]},
                            occupation=1.0, participation_ratio=0.1))
    cbm_info = EdgeInfo(band_idx=12, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=OrbitalInfo(
                            energy=1.2, orbitals={"Mn": [0.0, 0.0, 0.0, 0.0],
                                                  "O": [0.1, 0.3, 0.0, 0.0]},
                            occupation=0.0, participation_ratio=0.1))
    localized_orb_1 = LocalizedOrbital(
        band_idx=10, ave_energy=0.0, occupation=1.0,
        orbitals={"Mn": [0.1, 0.2, 0.0, 0.0], "O": [0.3, 0.4, 0.0, 0.0]},
        participation_ratio=0.1)
    localized_orb_2 = LocalizedOrbital(
        band_idx=11, ave_energy=1.0, occupation=0.0,
        orbitals={"Mn": [0.5, 0.8, 0.0, 0.0], "O": [0.0, 0.0, 0.0, 0.0]},
        participation_ratio=0.1)
    expected = BandEdgeStates(
        states=[BandEdgeState(vbm_info=vbm_info,
                              cbm_info=cbm_info,
                              vbm_orbital_diff=0.09999999999999998,
                              cbm_orbital_diff=0.09999999999999998,
                              localized_orbitals=[localized_orb_1,
                                                  localized_orb_2])])
    assert actual == expected


def test_make_band_edge_state_w_defect_charge_info(p_edge_state, orb_infos, mocker):
    defect_charge_info = mocker.Mock()
    defect_charge_info.localized_orbitals.return_value = [[9]]
    actual = make_band_edge_states(orb_infos, p_edge_state, defect_charge_info)
    expected = EdgeInfo(band_idx=8, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=OrbitalInfo(
                            energy=-1.1, orbitals={"Mn": [0.5, 0.6, 0.0, 0.0],
                                                   "O": [0.0, 0.0, 0.0, 0.0]},
                            occupation=1.0, participation_ratio=0.1))
    assert actual.states[0].vbm_info == expected


def test_orbital_diff(mocker):
    orb_1 = {"Mn": [0.1, 0.0, 0.0, 0.0]}
    orb_2 = {"Mn": [0.0, 0.1, 0.0, 0.0]}
    assert orbital_diff(orb_1, orb_2) == 0.2

    orb_1 = {"Mn": [0.1, 0.0, 0.0, 0.0]}
    orb_2 = {"O": [0.0, 0.1, 0.0, 0.0]}
    assert orbital_diff(orb_1, orb_2) == 0.2


