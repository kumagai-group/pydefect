# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest

from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues, \
    BandEdgeStates, EdgeState, EdgeCharacter, EdgeCharacters
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip


@pytest.fixture
def band_edge_eigenvalues():
    return BandEdgeEigenvalues(energies=np.array([[[0.0]]]),
                               occupations=np.array([[[0.0]]]),
                               lowest_band_energy=0.0,
                               highest_band_energy=1.0,
                               kpt_coords=[(0.0, 0.0, 0.0)],
                               kpt_weights=[1.0])


def test_band_edge_eigenvalues_msonable(band_edge_eigenvalues):
    assert_msonable(band_edge_eigenvalues)


def test_band_edge_eigenvalues_to_json_file(band_edge_eigenvalues, tmpdir):
    assert_json_roundtrip(band_edge_eigenvalues, tmpdir)


@pytest.fixture
def band_edge_states():
    return BandEdgeStates([EdgeState.donor_phs, EdgeState.no_in_gap],
                          method="manual")


def test_band_edge_states_is_shallow(band_edge_states):
    assert band_edge_states.is_shallow is True
    actual = BandEdgeStates([EdgeState.localized_state, EdgeState.no_in_gap],
                       method="manual")
    assert actual.is_shallow is False


def test_band_edge_states_msonable(band_edge_states):
    assert_msonable(band_edge_states)


def test_band_edge_states_to_json_file(band_edge_states, tmpdir):
    assert_json_roundtrip(band_edge_states, tmpdir)


@pytest.fixture
def edge_characters():
    edge_character = EdgeCharacter(hob_p_ratio=0.1,
                        lub_p_ratio=0.2,
                        hob_energy=0.3,
                        lub_energy=0.4,
                        hob_orbitals={"Mn": [0.5, 0.6]},
                        lub_orbitals={"Mn": [0.7, 0.8]})
    return EdgeCharacters(edge_characters=[edge_character])


def test_edge_characters_msonable(edge_characters):
    assert_msonable(edge_characters)


def test_edge_characters_to_json_file(edge_characters, tmpdir):
    assert_json_roundtrip(edge_characters, tmpdir)


"""
TODO
-

DONE
"""