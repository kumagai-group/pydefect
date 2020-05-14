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


def test_supercell_info(site, monoclinic):
    supercell_info = SupercellInfo(monoclinic,
                                   [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                   sites={"H1": site})
    assert_msonable(supercell_info)


def test_supercell_info_distances(monoclinic_supercell):
    sites = {"H1": Site(wyckoff_letter="a", site_symmetry="mmm",
                        equivalent_atoms=[0, 1, 2, 3]),
             "He1": Site(wyckoff_letter="b", site_symmetry="mmm",
                         equivalent_atoms=[4, 5, 6, 7])}
    supercell_info = SupercellInfo(monoclinic_supercell,
                                   [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                   sites=sites)
    coords = supercell_info.coords
    assert coords["H1"] == {"H": [3.91], "He": [2.5, 3.0, 3.5]}
    assert coords["He1"] == {"H": [2.5, 3.0, 3.5], "He": [3.91]}


"""
TODO
- Evaluate coords at Site
- Implement print info
- Generate default defect.in

DONE
"""