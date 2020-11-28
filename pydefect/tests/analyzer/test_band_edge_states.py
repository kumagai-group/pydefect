# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues, \
    BandEdgeStates, EdgeState, EdgeCharacter, EdgeCharacters
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip


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
def band_edge_states():
    return BandEdgeStates([EdgeState.donor_phs, EdgeState.no_in_gap])


def test_band_edge_states_is_shallow(band_edge_states):
    assert band_edge_states.is_shallow is True
    actual = BandEdgeStates([EdgeState.in_gap_state, EdgeState.no_in_gap])
    assert actual.is_shallow is False


def test_band_edge_states_msonable(band_edge_states):
    assert_msonable(band_edge_states)


def test_band_edge_states_to_json_file(band_edge_states, tmpdir):
    assert_json_roundtrip(band_edge_states, tmpdir)


def test_band_edge_states_yaml(band_edge_states, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    band_edge_states.to_yaml()
    expected = """up: Donor PHS
down: No in-gap state"""
    assert Path("band_edge_states.yaml").read_text() == expected

    actual = band_edge_states.from_yaml("band_edge_states.yaml")
    assert actual == band_edge_states


@pytest.fixture
def edge_characters():
    edge_character = EdgeCharacter(vbm_participation_ratio=0.1,
                                   cbm_participation_ratio=0.2,
                                   hob_bottom_e=0.3,
                                   lub_top_e=0.4,
                                   vbm=0.5,
                                   cbm=0.6,
                                   vbm_orbitals={"Mn": [0.5, 0.6]},
                                   cbm_orbitals={"Mn": [0.7, 0.8]})
    return EdgeCharacters([edge_character])


def test_edge_characters_msonable(edge_characters):
    assert_msonable(edge_characters)


def test_edge_characters_to_json_file(edge_characters, tmpdir):
    assert_json_roundtrip(edge_characters, tmpdir)


def test_str(edge_characters):
    assert str(edge_characters) == """  up: EdgeCharacter(hob_bottom_e=0.3, lub_top_e=0.4, vbm=0.5, cbm=0.6, vbm_orbitals={'Mn': [0.5, 0.6]}, cbm_orbitals={'Mn': [0.7, 0.8]}, vbm_participation_ratio=0.1, cbm_participation_ratio=0.2)"""


