# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues, \
    BandEdgeStates, OrbitalInfo, BandEdgeOrbitalInfos, PerfectBandEdgeState, \
    EdgeInfo, BandEdgeState, LocalizedOrbital, pretty_orbital
from pydefect.defaults import defaults
from vise.tests.helpers.assertion import assert_msonable, assert_json_roundtrip


@pytest.fixture
def band_edge_eigenvalues():
    return BandEdgeEigenvalues(energies_and_occupations=[[[[0.0, 1.0]]]],
                               kpt_coords=[[0.0, 0.0, 0.0]],
                               lowest_band_index=10)


def test_band_edge_eigenvalues_msonable(band_edge_eigenvalues):
    assert_msonable(band_edge_eigenvalues)


def test_band_edge_eigenvalues_to_json_file(band_edge_eigenvalues, tmpdir):
    assert_json_roundtrip(band_edge_eigenvalues, tmpdir)


@pytest.fixture
def orbital_info():
    return OrbitalInfo(energy=1.0,
                       orbitals={"Mn": [0.5, 0.4, 0.0, 0.0]},
                       occupation=1.0,
                       participation_ratio=0.1)


def test_orbital_info_msonable(orbital_info):
    assert_msonable(orbital_info)


def test_pretty_orbital(orbital_info):
    assert pretty_orbital(orbital_info.orbitals) == "Mn-s: 0.50, Mn-p: 0.40"


@pytest.fixture
def band_edge_orbital_infos(orbital_info):
    return BandEdgeOrbitalInfos(orbital_infos=[[[orbital_info]]],
                                kpt_coords=[(0.0, 0.0, 0.0)],
                                kpt_weights=[1.0],
                                lowest_band_index=10,
                                fermi_level=0.5)


def test_kpt_idx(band_edge_orbital_infos):
    assert band_edge_orbital_infos.kpt_idx([0.0, 0.0, 0.0]) == 0


def test_band_edge_orbital_info_json_roundtrip(band_edge_orbital_infos, tmpdir):
    assert_json_roundtrip(band_edge_orbital_infos, tmpdir)


def test_band_edge_orbital_info_energies_occupations(band_edge_orbital_infos):
    actual = band_edge_orbital_infos.energies_and_occupations
    expected = [[[[1.0, 1.0]]]]
    assert actual == expected


def test_band_edge_orbital_info_repr():
    orbital_info_1 = OrbitalInfo(energy=1.23456789,
                                 orbitals={"Mn": [0.5, 0.4, 0.0, 0.0],
                                           "Fe": [0.0, 0.0, 0.0, 0.0]},
                                 occupation=0.99998,
                                 participation_ratio=0.1111111)
    orbital_info_2 = OrbitalInfo(energy=2.34567890,
                                 orbitals={"Mn": [0.1, 0.1, 0.0, 0.0],
                                           "Fe": [0.0, 0.0, 0.0, 0.0]},
                                 occupation=0.0001,
                                 participation_ratio=0.2222222)

    band_edge_orbital_infos =  \
        BandEdgeOrbitalInfos(orbital_infos=[[[orbital_info_1, orbital_info_2]]],
                             kpt_coords=[(0.0, 0.0, 0.0)],
                             kpt_weights=[1.0],
                             lowest_band_index=10,
                             fermi_level=0.5)

    actual = band_edge_orbital_infos.__str__()
    print(actual)
    expected = """ -- band-edge orbitals info
K-points info
Index  Coords                    Weight
1      ( 0.000,  0.000,  0.000)  1.000

Band info near band edges
Index  Kpoint index  Energy  Occupation  P-ratio  Orbital
11     1             1.23    1.0         0.1      Mn-s: 0.50, Mn-p: 0.40
--
12     1             2.35    0.0         0.2
--
"""
    assert actual == expected


@pytest.fixture
def perfect_band_edge_state(orbital_info):
    vbm_info = EdgeInfo(band_idx=10, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=orbital_info)
    return PerfectBandEdgeState(vbm_info=vbm_info, cbm_info=vbm_info)


def test_perfect_edge_state_json_roundtrip(perfect_band_edge_state, tmpdir):
    assert_json_roundtrip(perfect_band_edge_state, tmpdir)


@pytest.fixture
def band_edge_states(orbital_info):
    vbm_info = EdgeInfo(band_idx=10, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=orbital_info)
    cbm_info = EdgeInfo(band_idx=12, kpt_coord=(0.0, 0.0, 0.0),
                        orbital_info=orbital_info)
    localized_orbital = LocalizedOrbital(
        band_idx=11, ave_energy=1.9, occupation=1.0,
        orbitals={"H": [1.0, 0.0, 0.0, 0.0]})
    band_edge_state = BandEdgeState(vbm_info=vbm_info, cbm_info=cbm_info,
                                    vbm_orbital_diff=0.5,
                                    cbm_orbital_diff=0.5,
                                    localized_orbitals=[localized_orbital],
                                    vbm_hole_occupation=defaults.state_occupied_threshold - 1e-5,
                                    cbm_electron_occupation=defaults.state_occupied_threshold + 1e-5,
                                    )
    return BandEdgeStates(states=[band_edge_state])


def test_band_edge_states_json_roundtrip(band_edge_states, tmpdir):
    assert_json_roundtrip(band_edge_states, tmpdir)


def test_band_edge_states_band_indices(band_edge_states):
    assert band_edge_states.band_indices_from_vbm_to_cbm == [10, 11, 12]


def test_band_edge_states_has_donor_phs(band_edge_states):
    assert band_edge_states.has_donor_phs is True
    assert band_edge_states.has_acceptor_phs is False
    assert band_edge_states.has_unoccupied_localized_state is False
    assert band_edge_states.has_occupied_localized_state is True


def test_band_edge_states_str(band_edge_states):
    expected = f""" -- band-edge states info
Spin-up
     Index  Energy  P-ratio  Occupation  OrbDiff  Orbitals                K-point coords
VBM  11     1.000   0.10     1.00        0.50     Mn-s: 0.50, Mn-p: 0.40  ( 0.000,  0.000,  0.000)
CBM  13     1.000   0.10     1.00        0.50     Mn-s: 0.50, Mn-p: 0.40  ( 0.000,  0.000,  0.000)
vbm has acceptor phs: False ({defaults.state_occupied_threshold:.3f} vs. {defaults.state_occupied_threshold})
cbm has donor phs: True ({defaults.state_occupied_threshold:.3f} vs. {defaults.state_occupied_threshold})
---
Localized Orbital(s)
Index  Energy  P-ratio  Occupation  Orbitals
12     1.900   None     1.00        H-s: 1.00
"""
    print(band_edge_states)
    assert band_edge_states.__str__() == expected


def test_band_edge_states_str2(band_edge_states):
    copy = deepcopy(band_edge_states)
    copy.states[0].localized_orbitals = [LocalizedOrbital(
        band_idx=11, ave_energy=1.9, occupation=1.0,
        orbitals={"H": [1.0, 0.0, 0.0, 0.0]}, radius=1.000001,
        center=[0.1234, 0.0, 0.0])]
    actual = copy.__str__()
    print(actual)
    expected = f""" -- band-edge states info
Spin-up
     Index  Energy  P-ratio  Occupation  OrbDiff  Orbitals                K-point coords
VBM  11     1.000   0.10     1.00        0.50     Mn-s: 0.50, Mn-p: 0.40  ( 0.000,  0.000,  0.000)
CBM  13     1.000   0.10     1.00        0.50     Mn-s: 0.50, Mn-p: 0.40  ( 0.000,  0.000,  0.000)
vbm has acceptor phs: False ({defaults.state_occupied_threshold:.3f} vs. {defaults.state_occupied_threshold})
cbm has donor phs: True ({defaults.state_occupied_threshold:.3f} vs. {defaults.state_occupied_threshold})
---
Localized Orbital(s)
Index  Energy  P-ratio  Occupation  Orbitals   Radius  Center
12     1.900   None     1.00        H-s: 1.00  1.00    ( 0.123,  0.000,  0.000)
"""
    assert actual == expected