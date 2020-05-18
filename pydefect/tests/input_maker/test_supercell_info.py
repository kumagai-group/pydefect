# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pydefect.input_maker.supercell_info import Site
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def site():
    return Site(element="H",
                wyckoff_letter="a",
                site_symmetry="m3m",
                equivalent_atoms=[0, 1, 2, 3])


def test_msonable(site):
    assert_msonable(site)


def test_supercell_info(supercell_info):
    assert_msonable(supercell_info)


def test_supercell_info_distances(supercell_info):
    assert supercell_info.coords("H1") == ({"H": [3.91], "He": [2.5, 3.0, 3.5]}, 4.25)
    assert supercell_info.coords("He1") == ({"H": [2.5, 3.0, 3.5], "He": [3.91]}, 4.25)


def test_supercell_info_str(supercell_info):
    expected = """  Space group: Fmmm
Transformation matrix: [1, 0, 0]  [0, 1, 0]  [0, 0, 1]
Cell multiplicity: 1

   Irreducible element: H1
        Wyckoff letter: a
         Site symmetry: mmm
         Cutoff radius: 4.25
          Coordination: {'H': [3.91], 'He': [2.5, 3.0, 3.5]}
      Equivalent atoms: [0, 1, 2, 3]
Fractional coordinates: 0.0000000  0.0000000  0.0000000
     Electronegativity: 2.2
       Oxidation state: 1

   Irreducible element: He1
        Wyckoff letter: b
         Site symmetry: mmm
         Cutoff radius: 4.25
          Coordination: {'H': [2.5, 3.0, 3.5], 'He': [3.91]}
      Equivalent atoms: [4, 5, 6, 7]
Fractional coordinates: 0.0000000  0.0000000  0.5000000
     Electronegativity: None
       Oxidation state: 0
"""
    assert str(supercell_info) == expected


def test_from_json(tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    actual = supercell_info.from_json_file().as_dict()
    assert actual == supercell_info.as_dict()


"""
TODO
- Implement from_json

DONE
- Evaluate coords at Site
- Add spacegroup to Supercell Info
- Implement print info
"""