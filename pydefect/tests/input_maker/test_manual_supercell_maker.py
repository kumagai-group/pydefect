# -*- coding: utf-8 -*-
#  Copyright (c) 2023. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from pydefect.input_maker.manual_supercell_maker import ManualSupercellMaker, \
    make_sites_from_yaml_file, IncorrectSiteYamlFileError, \
    check_site_structure_consistency
from pydefect.input_maker.supercell_info import SimpleSite, SupercellInfo


def test_create_supercell(simple_cubic_2x1x1, simple_cubic_2x2x2):
    matrix = [[1, 0, 0], [0, 2, 0], [0, 0, 2]]
    sites = {"H1": SimpleSite("H", 1)}
    actual = ManualSupercellMaker(simple_cubic_2x1x1,
                                  matrix=matrix, sites=sites)
    assert actual.supercell_info.structure == simple_cubic_2x2x2
    assert actual.supercell_info.sites == {"H1": SimpleSite("H", 4)}


def test_create_supercell_2(simple_cubic, simple_cubic_2x2x2):
    matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    sites = {"H1": SimpleSite("H", 0)}
    actual = ManualSupercellMaker(simple_cubic,
                                  sites=sites,
                                  min_num_atoms=5,
                                  max_num_atoms=10)
    expected = SupercellInfo(simple_cubic_2x2x2,
                             "Pm-3m",
                             matrix,
                             sites,
                             unitcell_structure=simple_cubic)
    assert actual.supercell_info == expected


def test_make_sites_from_yaml(tmpdir):
    tmpdir.chdir()
    Path("sites.yaml").write_text("""
H1:
  site_index: 1
  site_symmetry: P1
""")
    actual = make_sites_from_yaml_file("sites.yaml")
    expected = {"H1": SimpleSite("H", 1, "P1")}
    assert actual == expected


def test_make_sites_from_yaml_w_error(tmpdir):
    tmpdir.chdir()
    Path("sites.yaml").write_text("""
H1:
  atom: H
""")
    with pytest.raises(TypeError):
        make_sites_from_yaml_file("sites.yaml")


def test_check_site_structure_consistency(simple_cubic):
    sites = {"H": SimpleSite("H", 0)}
    assert check_site_structure_consistency(sites, simple_cubic) is None

    sites = {"He": SimpleSite("He", 0)}
    with pytest.raises(IncorrectSiteYamlFileError):
        check_site_structure_consistency(sites, simple_cubic)

    sites = {"H": SimpleSite("H", 1)}
    with pytest.raises(IncorrectSiteYamlFileError):
        check_site_structure_consistency(sites, simple_cubic)
