# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import EdgeCharacter, EdgeState
from pydefect.analyzer.make_band_edge_state import make_band_edge_state, \
    are_orbitals_similar


@pytest.fixture
def no_in_gap_character():
    return EdgeCharacter(vbm_participation_ratio=0.1,
                         cbm_participation_ratio=0.2,
                         hob_bottom_e=-0.1,
                         lub_top_e=1.1,
                         vbm=0.0,
                         cbm=1.0,
                         vbm_orbitals={"Mn": [0.5, 0.6]},
                         cbm_orbitals={"O": [0.7, 0.8]})


@pytest.fixture
def donor_phs_character():
    return EdgeCharacter(vbm_participation_ratio=0.1,
                         cbm_participation_ratio=0.2,
                         hob_bottom_e=0.9,
                         lub_top_e=1.1,
                         vbm=None,
                         cbm=None,
                         vbm_orbitals={"O": [0.5, 0.6]},
                         cbm_orbitals={"O": [0.7, 0.8]})


@pytest.fixture
def acceptor_phs_character():
    return EdgeCharacter(vbm_participation_ratio=0.1,
                         cbm_participation_ratio=0.2,
                         hob_bottom_e=-0.2,
                         lub_top_e=0.1,
                         vbm=None,
                         cbm=None,
                         vbm_orbitals={"Mn": [0.5, 0.6]},
                         cbm_orbitals={"Mn": [0.7, 0.8]})


@pytest.fixture
def in_gap_state_character():
    return EdgeCharacter(vbm_participation_ratio=0.5,
                         cbm_participation_ratio=0.5,
                         hob_bottom_e=0.4,
                         lub_top_e=0.6,
                         vbm=0.45,
                         cbm=0.55,
                         vbm_orbitals={"Mn": [0.3, 0.2], "O": [0.3, 0.1]},
                         cbm_orbitals={"Mn": [0.3, 0.2], "O": [0.4, 0.2]})


# @pytest.fixture
# def unknown_state_character():
#     return EdgeCharacter(vbm_participation_ratio=0.1,
#                          cbm_participation_ratio=0.2,
#                          hob_bottom_e=0.4,
#                          lub_top_e=0.6,
#                          vbm=0.45,
#                          cbm=0.55,
#                          vbm_orbitals={"Mn": [0.3, 0.2], "O": [0.3, 0.1]},
#                          cbm_orbitals={"Mn": [0.3, 0.2], "O": [0.4, 0.2]})


@pytest.fixture
def perfect_character():
    return EdgeCharacter(hob_bottom_e=-0.1,
                         lub_top_e=1.1,
                         vbm=0.0,
                         cbm=1.0,
                         vbm_orbitals={"Mn": [0.5, 0.6]},
                         cbm_orbitals={"O": [0.7, 0.8]})


def test_make_band_edge_state_no_in_gap(no_in_gap_character, perfect_character):
    actual = make_band_edge_state(no_in_gap_character, perfect_character)
    assert actual == EdgeState.no_in_gap


def test_make_band_edge_state_donor_phs(donor_phs_character, perfect_character):
    actual = make_band_edge_state(donor_phs_character, perfect_character)
    assert actual == EdgeState.donor_phs


def test_make_band_edge_state_acceptor_phs(acceptor_phs_character,
                                           perfect_character):
    actual = make_band_edge_state(acceptor_phs_character, perfect_character)
    assert actual == EdgeState.acceptor_phs


def test_make_band_edge_state_in_gap_state(in_gap_state_character,
                                           perfect_character):
    actual = make_band_edge_state(in_gap_state_character, perfect_character)
    assert actual == EdgeState.in_gap_state


# def test_make_band_edge_state_unknown_state(unknown_state_character, perfect_character):
#     actual = make_band_edge_state(unknown_state_character, perfect_character)
#     assert actual == EdgeState.unknown


def test_are_orbitals_similar(mocker):
    mock = mocker.patch("pydefect.analyzer.make_band_edge_state.defaults")
    mock.similar_orb_criterion = 27.0  # = 1 + 2 + 6 + 7 + 8 + 1 + 1 + 1
    orb_1 = {"Mn": [0.0, 1.0, 2.0], "O": [3.0, 4.0, 5.0]}
    orb_2 = {"Mg": [-6.0, -7.0, -8.0], "O": [4.0, 3.0, 6.0]}
    assert are_orbitals_similar(orb_1, orb_2) is False

    mock.similar_orb_criterion = 27.1
    assert are_orbitals_similar(orb_1, orb_2) is True


dir_path = Path(__file__).parent / "band_edge_database" / "NaCl"


def test_edge_character_nacl_acceptor_phs():
    ref = loadfn(dir_path / "perfect_edge_characters.json")[0]  # only spin up
    target = loadfn(dir_path / "Va_Na1_1_edge_characters.json")[0]
    assert make_band_edge_state(target, ref) == EdgeState.acceptor_phs


def test_edge_character_nacl_no_in_gap():
    ref = loadfn(dir_path / "perfect_edge_characters.json")[0]
    target = loadfn(dir_path / "Va_Na1_-1_edge_characters.json")[0]
    assert make_band_edge_state(target, ref) == EdgeState.no_in_gap


