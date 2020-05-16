# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pymatgen import Structure, Lattice
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.defaults import defaults

# "H" at [0.0, 0.0, 0.0] is removed here.
coords = \
    [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
     [0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.5]]
perturbed_coords = \
    [[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
     [0.0, 0.0, 0.501], [0.0, 0.501, 0.0], [0.501, 0.0, 0.0], [0.5, 0.5, 0.5]]

rocksalt = Structure(
    Lattice.cubic(10.0), ["H"] * 3 + ["He"] * 4, coords)
perturbed_rocksalt = Structure(
    Lattice.cubic(10.0), ["H"] * 3 + ["He"] * 4, perturbed_coords)


@pytest.fixture
def defect_entry():
    return DefectEntry(name="Va_O1",
                       charge=1,
                       initial_structure=rocksalt,
                       perturbed_initial_structure=perturbed_rocksalt,
                       initial_site_symmetry="m-3m")


def test_defect_entry(defect_entry):
    assert defect_entry.name == "Va_O1"




"""
TODO
- Store basic info for a single input defect initial structure.
- create prior_info
- from structure
o                name: str,
?                defect_type: DefectType,
o                initial_structure: Structure,
o                perturbed_initial_structure: Structure,
o                removed_atoms: list,
o                inserted_atoms: list,
x                changes_of_num_elements: dict,
o                charge: int,
o                initial_site_symmetry: str,
o                cutoff: float,
x                neighboring_sites: list,
x                annotation: Optional[str] = None,
x                multiplicity: Optional[int] = None):
DONE
"""