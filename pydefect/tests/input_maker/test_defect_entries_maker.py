# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np

from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_entries_maker import (
    DefectEntriesMaker, copy_to_structure, to_istructure, add_atom_to_structure,
    perturb_structure, random_3d_vector)
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.defect_set import DefectSet
from pydefect.defaults import defaults
from pymatgen import Site


def test_defect_entries_maker(supercell_info, ortho_conventional, mocker):
    mock = mocker.patch("pydefect.input_maker.defect_entries_maker.defaults")
    mock.displace_distance = 0.0
    defect_set = DefectSet(defects={SimpleDefect(None, "H1", [0, 1]),
                                    SimpleDefect("Ne", "He1", [0])})
    maker = DefectEntriesMaker(supercell_info, defect_set)

    va_h1_str = copy_to_structure(ortho_conventional)
    va_h1_str.pop(0)
    va_h1_str = to_istructure(va_h1_str)

    ne_he1_str = copy_to_structure(ortho_conventional)
    coords = ne_he1_str.pop(4).frac_coords
    ne_he1_str.append("Ne", coords)
    ne_he1_str = to_istructure(ne_he1_str)

    defect_entries = {
        DefectEntry("Ne_He1", 0, ne_he1_str, ne_he1_str.copy(), "mmm"),
        DefectEntry("Va_H1", 0, va_h1_str, va_h1_str.copy(), "mmm"),
        DefectEntry("Va_H1", 1, va_h1_str, va_h1_str.copy(), "mmm"),
    }

    assert maker.defect_entries == defect_entries


def test_defect_entries_maker_insert_host_atoms(supercell_info, mocker):
    mock = mocker.patch("pydefect.input_maker.defect_entries_maker.defaults")
    mock.displace_distance = 0.0

    defect_set = DefectSet(defects={SimpleDefect("He", "H1", [0]),
                                    SimpleDefect("H", "He1", [0])})
    maker = DefectEntriesMaker(supercell_info, defect_set)

    for defect_entry in maker.defect_entries:
        if defect_entry.name == "H_He1":
            actual = defect_entry.initial_structure[0].frac_coords
            np.testing.assert_array_equal(actual, [0.0, 0.0, 0.5])
        if defect_entry.name == "He_H1":
            print(defect_entry.initial_structure)
            actual = defect_entry.initial_structure[3].frac_coords
            np.testing.assert_array_equal(actual, [0.0, 0.0, 0.0])


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
    mock = mocker.patch("pydefect.input_maker.defect_entries_maker.random_3d_vector")
    mock.return_value = [0.0, 0.0, 0.35]
    structure = copy_to_structure(ortho_conventional)
    actual = perturb_structure(structure,
                               center=[0.0, 0.0, 0.0], cutoff=2.5-1e-5)
    expected = structure
    expected[0].frac_coords = [0.0, 0.0, 0.05]

    mock.assert_called_with(defaults.displace_distance)
    assert actual == expected


def test_random_3d_vector():
    actual = random_3d_vector(1.0)
    assert np.linalg.norm(actual) <= 1.0
    assert len(actual) == 3

"""
TODO
- Take SupercellInfo and DefectSet, return List[DefectEntry]
- Insert host antisite atoms to the middle of structure 
- Set cutoff  properly.

DONE
"""