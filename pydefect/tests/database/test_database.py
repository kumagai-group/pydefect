# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.database.database import electronegativity, oxidation_state, rcore


def test_electronegativity():
    assert electronegativity("H") == 2.20


def test_electronegativity_not_exist():
    assert electronegativity("He") is None


def test_oxidation_state():
    assert oxidation_state("H") == 1


def test_oxidation_state_not_exist():
    assert oxidation_state("Os") == 0


def test_rcore():
    assert rcore["H"] == 0.582
