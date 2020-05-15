# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen import Element
from pydefect.input_maker.supercell_info import Site, SupercellInfo
from pydefect.tests.helpers.assertion import assert_msonable
from vise.util.structure_symmetrizer import StructureSymmetrizer


@pytest.fixture
def site():
    return Site(wyckoff_letter="a",
                site_symmetry="m3m",
                equivalent_atoms=[0, 1, 2, 3])


def test_site(site):
    assert_msonable(site)


@pytest.fixture
def supercell_info(ortho_conventional):
    sites = {"H1": Site(wyckoff_letter="a", site_symmetry="mmm",
                        equivalent_atoms=[0, 1, 2, 3]),
             "He1": Site(wyckoff_letter="b", site_symmetry="mmm",
                         equivalent_atoms=[4, 5, 6, 7])}
    return SupercellInfo(ortho_conventional,
                         "Fmmm",
                         [[1, 0, 0], [0, 1, 0], [0, 0, 1]], sites)
#
# def test(ortho_conventional):
#     print(StructureSymmetrizer(ortho_conventional).spglib_sym_data)


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
     Electronegativity: 0.0
       Oxidation state: 0
"""
    assert str(supercell_info) == expected


"""
TODO
- Add spacegroup to Supercell Info
- Implement print info
- Generate default defect.in

DONE
- Evaluate coords at Site
"""