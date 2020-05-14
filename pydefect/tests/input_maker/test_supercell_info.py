# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pymatgen import Element
from pydefect.input_maker.supercell_info import Site, SupercellInfo
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def site():
    return Site(element=Element.H,
                wyckoff_letter="a",
                site_symmetry="m3m",
                cutoff_radius=2.0,
                coordination={"O": [1.2, 1.3]},
                equivalent_atoms=[0, 1, 2, 3])


def test_site(site):
    assert_msonable(site)


def test_supercell_info(site, monoclinic):
    supercell_info = SupercellInfo(monoclinic,
                                   [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                   sites={"H1": site})
    assert_msonable(supercell_info)


"""
TODO
- Hold some info

DONE
"""