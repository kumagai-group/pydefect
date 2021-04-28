# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from copy import deepcopy

import pytest
from pydefect.input_maker.local_extrema import VolumetricDataLocalExtrema, \
    CoordInfo, VolumetricDataAnalyzeParams
from pydefect.input_maker.supercell_info import Interstitial, SupercellInfo
from pydefect.util.structure_tools import Coordination
from pymatgen import Structure
from vise.tests.helpers.assertion import assert_json_roundtrip
from vise.util.structure_symmetrizer import Site


@pytest.fixture
def local_extrema(simple_cubic, vol_params):
    coordination = Coordination({"Mn": [1.0, 2.0]}, cutoff=4.0,
                                neighboring_atom_indices=[1, 2])
    local_extremum = CoordInfo(site_symmetry="3m",
                               coordination=coordination,
                               frac_coords=[(0.1, 0.1, 0.1)],
                               quantities=[1.234])

    coordination = Coordination({"Mn": [1.0, 2.0]}, cutoff=4.0,
                                neighboring_atom_indices=[1, 2])
    local_extremum2 = CoordInfo(site_symmetry="4mm",
                                coordination=coordination,
                                frac_coords=[(0.1, 0.0, 0.0)],
                                quantities=[2.101])

    return VolumetricDataLocalExtrema(unit_cell=simple_cubic, is_min=True,
                                      extrema_points=[local_extremum,
                                                      local_extremum2],
                                      info="test",
                                      params=vol_params)


def test_local_extrema_json_roundtrip(local_extrema, tmpdir):
    tmpdir.chdir()
    assert_json_roundtrip(local_extrema, tmpdir)


def test_local_extrema_str(local_extrema, tmpdir):
    actual = local_extrema.__str__()
    expected = """info: test
min_or_max: min
extrema_points:
#  site_sym  coordination        frac_coords               quantity
1  3m        {'Mn': [1.0, 2.0]}  ( 0.100,  0.100,  0.100)  1.2
2  4mm       {'Mn': [1.0, 2.0]}  ( 0.100,  0.000,  0.000)  2.1"""
    print(actual)
    assert actual == expected


def test_local_extrema_append_sites_to_supercell_info(local_extrema,
                                                      simple_cubic):
    sites = {"H1": Site(element="H", wyckoff_letter="a", site_symmetry="mmm",
                        equivalent_atoms=[0, 1, 2, 3, 4, 5, 6, 7])}
    actual_cubic = Structure.from_dict(simple_cubic.as_dict())
    supercell_info = SupercellInfo(actual_cubic * [2, 2, 2],
                                   "Pm-3m",
                                   [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                                   sites, [])
    actual = local_extrema.append_sites_to_supercell_info(supercell_info, indices=[1, 2])

    interstitials = [Interstitial(frac_coords=[0.05]*3, site_symmetry="3m", info="test #1"),
                     Interstitial(frac_coords=[0.05, 0.0, 0.0], site_symmetry="4mm", info="test #2")]
    expected = SupercellInfo(actual_cubic * [2, 2, 2],
                             "Pm-3m",
                             [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                             sites, interstitials)
    assert actual == expected


def test_local_extrema_append_sites_to_supercell_info2(local_extrema,
                                                       simple_cubic):
    sites = {"H1": Site(element="H", wyckoff_letter="a", site_symmetry="mmm",
                        equivalent_atoms=[0, 1, 2, 3, 4, 5, 6, 7])}
    actual_cubic = Structure.from_dict(simple_cubic.as_dict())
    supercell_info = SupercellInfo(actual_cubic * [2, 2, 2],
                                   "Pm-3m",
                                   [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                                   sites, [])

    actual_2 = local_extrema.append_sites_to_supercell_info(supercell_info, indices=[2])

    interstitials = [Interstitial(frac_coords=[0.05, 0.0, 0.0], site_symmetry="4mm", info="test #2")]
    expected_2 = SupercellInfo(actual_cubic * [2, 2, 2],
                             "Pm-3m",
                             [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                             sites, interstitials)
    assert actual_2 == expected_2
