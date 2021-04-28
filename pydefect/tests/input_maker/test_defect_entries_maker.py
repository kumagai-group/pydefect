# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest

from pydefect.defaults import defaults
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entries_maker import (
    DefectEntriesMaker, copy_to_structure, to_istructure, add_atom_to_structure,
    perturb_structure, random_3d_vector)
from pydefect.input_maker.defect_entry import DefectEntry, PerturbedSite
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.supercell_info import SupercellInfo, Site, \
    Interstitial
from pymatgen import PeriodicSite
from vise.tests.helpers.assertion import assert_dataclass_almost_equal


@pytest.fixture
def cubic_supercell_info(cubic_supercell):
    sites = {"H1": Site(element="H", wyckoff_letter="a", site_symmetry="m-3m",
                        equivalent_atoms=list(range(32))),
             "He1": Site(element="He", wyckoff_letter="b", site_symmetry="m-3m",
                         equivalent_atoms=list(range(32, 64)))}
    interstitial = Interstitial([0.25]*3, site_symmetry="yy")
    return SupercellInfo(cubic_supercell,
                         "Fm-3m",
                         [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                         sites,
                         interstitials=[interstitial])


def test_defect_entries_maker(cubic_supercell_info, cubic_supercell, mocker):
    mock = mocker.patch("pydefect.input_maker.defect_entries_maker.defaults")
    mock.displace_distance = 0.0
    defect_set = DefectSet(defects={SimpleDefect(None, "H1", [1]),
                                    SimpleDefect("Ne", "He1", [0])})
    maker = DefectEntriesMaker(cubic_supercell_info, defect_set)

    va_h1_str = cubic_supercell.copy()
    va_h1_str.pop(0)
    va_h1_str = to_istructure(va_h1_str)

    ne_he1_str = cubic_supercell.copy()
    coords = ne_he1_str.pop(32).frac_coords
    ne_he1_str.append("Ne", coords)
    ne_he1_str = to_istructure(ne_he1_str)

    defect_entries = {DefectEntry(name="Ne_He1",
                                  charge=0,
                                  structure=ne_he1_str,
                                  site_symmetry="m-3m",
                                  defect_center=(0.25, 0.0, 0.0)),
                      DefectEntry(name="Va_H1",
                                  charge=1,
                                  structure=va_h1_str,
                                  site_symmetry="m-3m",
                                  defect_center=(0.0, 0.0, 0.0))}

    assert maker.defect_entries == defect_entries


def test_defect_entries_maker_insert_host_atoms(cubic_supercell_info, mocker):
    mock = mocker.patch("pydefect.input_maker.defect_entries_maker.defaults")
    mock.displace_distance = 0.0

    defect_set = DefectSet(defects={SimpleDefect("He", "H1", [0]),
                                    SimpleDefect("H", "He1", [0])})
    maker = DefectEntriesMaker(cubic_supercell_info, defect_set)

    for defect_entry in maker.defect_entries:
        if defect_entry.name == "H_He1":
            actual = defect_entry.structure[0].frac_coords
            np.testing.assert_array_equal(actual, [0.25, 0.0, 0.0])
        if defect_entry.name == "He_H1":
            actual = defect_entry.structure[31].frac_coords
            np.testing.assert_array_equal(actual, [0.0, 0.0, 0.0])


def test_defect_entries_interstitials(cubic_supercell_info, mocker):
    mock = mocker.patch("pydefect.input_maker.defect_entries_maker.defaults")
    mock.displace_distance = 0.0

    defect_set = DefectSet(defects={SimpleDefect("H", "i1", [0])})
    maker = DefectEntriesMaker(cubic_supercell_info, defect_set)

    defect_entry = maker.defect_entries.pop()
    actual = defect_entry.structure[0].frac_coords
    np.testing.assert_array_equal(actual, [0.25, 0.25, 0.25])


def test_add_atom_to_structure(ortho_conventional):
    coords = [0.1, 0.1, 0.1]

    structure_1 = copy_to_structure(ortho_conventional)
    add_atom_to_structure(structure_1, "H", coords)
    assert len(structure_1) == 9
    np.testing.assert_array_equal(structure_1[0].frac_coords, coords)

    structure_2 = copy_to_structure(ortho_conventional)
    add_atom_to_structure(structure_2, "He", coords)
    assert len(structure_2) == 9
    np.testing.assert_array_equal(structure_2[4].frac_coords, coords)

    structure_3 = copy_to_structure(ortho_conventional)
    add_atom_to_structure(structure_3, "Li", coords)
    assert len(structure_3) == 9
    np.testing.assert_array_equal(structure_3[8].frac_coords, coords)


def test_perturbed_structure(ortho_conventional, mocker):
    mock = mocker.patch(
        "pydefect.input_maker.defect_entries_maker.random_3d_vector")
    mock.return_value = ([0.0, 0.0, 0.35], 1.0)
    structure = copy_to_structure(ortho_conventional)
    actual = perturb_structure(structure,
                               center=[0.0, 0.0, 0.0], cutoff=2.5-1e-5)
    expected_structure = structure
    expected_structure[0].frac_coords = [0.0, 0.0, 0.05]
    expected_site = PerturbedSite(element='H',
                                  distance=0.0,
                                  initial_coords=(0.0, 0.0, 0.0),
                                  perturbed_coords=(0.0, 0.0, 0.05),
                                  displacement=1.0)

    mock.assert_called_with(defaults.displace_distance)
    assert actual[0] == expected_structure
    assert_dataclass_almost_equal(actual[1][0], expected_site)


def test_random_3d_vector():
    actual_vec, actual_dist = random_3d_vector(1.0)
    assert actual_dist <= 1.0
    assert len(actual_vec) == 3

    actual_vec, actual_dist = random_3d_vector(0.0)
    assert actual_dist == 0.0

"""
TODO
- Set cutoff  properly.

DONE
- Take SupercellInfo and DefectSet, return List[DefectEntry]
- Insert host antisite atoms to the middle of structure 
"""