# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues, \
    BandEdgeStates, OrbitalInfo, BandEdgeOrbitalInfos, PerfectBandEdgeState, \
    EdgeInfo, BandEdgeState, LocalizedOrbital, IsShallow
from vise.tests.helpers.assertion import assert_msonable, assert_json_roundtrip


@pytest.fixture
def band_edge_eigenvalues():
    return BandEdgeEigenvalues(energies_and_occupations=[[[[0.0, 1.0]]]],
                               kpt_coords=[(0.0, 0.0, 0.0)],
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


def test_pretty_orbitals(orbital_info):
    assert orbital_info.pretty_orbital() == "Mn-s: 0.50, Mn-p: 0.40"


@pytest.fixture
def band_edge_orbital_infos(orbital_info):
    return BandEdgeOrbitalInfos(orbital_infos=[[[orbital_info]]],
                                kpt_coords=[(0.0, 0.0, 0.0)],
                                kpt_weights=[1.0],
                                lowest_band_index=10,
                                fermi_level=0.5)


def test_band_edge_orbital_info_json_roundtrip(band_edge_orbital_infos, tmpdir):
    assert_json_roundtrip(band_edge_orbital_infos, tmpdir)


def test_band_edge_orbital_info_energies_occupations(band_edge_orbital_infos):
    actual = band_edge_orbital_infos.energies_and_occupations
    expected = [[[[1.0, 1.0]]]]
    assert actual == expected


def test_band_edge_orbital_info_repr(band_edge_orbital_infos):
    actual = band_edge_orbital_infos.__repr__()
    print(actual)
    expected = """k-points
---  ---------------  ------
idx  coords           weight
1    (0.0, 0.0, 0.0)  1.0
---  ---------------  ------
--------  -----  ------  -----  -------  ----------------------
band_idx  k_idx  energy  occup  p_ratio  orbital
11        1      1.0     1.0    0.1      Mn-s: 0.50, Mn-p: 0.40
--

--------  -----  ------  -----  -------  ----------------------"""
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
                                    localized_orbitals=[localized_orbital])
    return BandEdgeStates(states=[band_edge_state])


def test_band_edge_states_json_roundtrip(band_edge_states, tmpdir):
    assert_json_roundtrip(band_edge_states, tmpdir)


@pytest.fixture
def is_shallow():
    return IsShallow(is_shallow=True)


def test_is_shallow_yaml(tmpdir):
    tmpdir.chdir()
    IsShallow(True).to_yaml()
    assert IsShallow.from_yaml() == IsShallow(True)
    Path("a.yaml").write_text("t")
    assert IsShallow.from_yaml("a.yaml") == IsShallow(True)
    Path("a.yaml").write_text("f")
    assert IsShallow.from_yaml("a.yaml") == IsShallow(False)
