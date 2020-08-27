# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.input_maker.supercell_info import Interstitial, SupercellInfo
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip
from pydefect.util.structure_tools import Coordination


@pytest.fixture
def interstitial():
    return Interstitial(frac_coords=[0.25, 0.25, 0.25],
                        wyckoff_letter="a",
                        site_symmetry="m3m")


def test_interstitial(interstitial):
    assert_msonable(interstitial)


def test_supercell_info_msonable(supercell_info):
    assert_msonable(supercell_info)


def test_supercell_info_from_dict(supercell_info):
    d = supercell_info.as_dict()
    d["transform_matrix"] = d.pop("transformation_matrix")
    expected = SupercellInfo.from_dict(d)
    assert expected == supercell_info


def test_supercell_info_distances(supercell_info):
    assert supercell_info.coords("H1") == Coordination({"H": [3.91], "He": [2.5, 3.0, 3.5]}, 4.25, {1, 4, 5, 6})
    assert supercell_info.coords("He1") == Coordination({"H": [2.5, 3.0, 3.5], "He": [3.91]}, 4.25, {0, 2, 3, 7})


def test_supercell_info_interstitials_distances(supercell_info):
    assert supercell_info.interstitial_coords(0) == Coordination({'H': [2.62, 2.62, 2.62, 2.62], 'He': [2.62, 2.62, 2.62, 2.62]}, 4.457, {0, 1, 2, 3, 4, 5, 6, 7})


def test_supercell_info_str(supercell_info):
    expected = """Space group: Fmmm
Transformation matrix: [1, 0, 0]  [0, 1, 0]  [0, 0, 1]
Cell multiplicity: 1

   Irreducible element: H1
        Wyckoff letter: a
         Site symmetry: mmm
         Cutoff radius: 4.25
          Coordination: {'H': [3.91], 'He': [2.5, 3.0, 3.5]}
      Equivalent atoms: 0..3
Fractional coordinates: 0.0000000  0.0000000  0.0000000
     Electronegativity: 2.2
       Oxidation state: 1

   Irreducible element: He1
        Wyckoff letter: b
         Site symmetry: mmm
         Cutoff radius: 4.25
          Coordination: {'H': [2.5, 3.0, 3.5], 'He': [3.91]}
      Equivalent atoms: 4..7
Fractional coordinates: 0.0000000  0.0000000  0.5000000
     Electronegativity: None
       Oxidation state: 0

-- interstitials
#1
Fractional coordinates: 0.2500000  0.2500000  0.2500000
        Wyckoff letter: x
         Site symmetry: yy
          Coordination: {'H': [2.62, 2.62, 2.62, 2.62], 'He': [2.62, 2.62, 2.62, 2.62]}
"""
    assert str(supercell_info) == expected


def test_json_round_trip(supercell_info, tmpdir):
    assert_json_roundtrip(supercell_info, tmpdir)


