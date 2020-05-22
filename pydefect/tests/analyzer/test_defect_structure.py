# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from pymatgen import Structure

from pydefect.analyzer.defect_structure import symmetrize_defect_structure, \
    DefectStructure
from pydefect.tests.helpers.assertion import assert_msonable


@pytest.fixture
def defect_structure(cubic_supercell_w_vacancy,
                     cubic_supercell_w_vacancy_w_perturb):
    return DefectStructure(name="test",
                           charge=1,
                           initial_structure=cubic_supercell_w_vacancy,
                           final_structure=cubic_supercell_w_vacancy_w_perturb,
                           initial_symmetry="m-3m",
                           final_symmetry="4mm",
                           defect_center=[0.0, 0.0, 0.0])


def test_defect_structure_msonable(defect_structure):
    assert_msonable(defect_structure)


def test_defect_structure_anchor_atom_index(defect_structure):
    assert defect_structure.anchor_atom_index == 6


def test_symmetrize_defect_structure():
    structure = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
5 0 0
0 5 0
0 0 5
Mg   O
4     3
Direct
0.0051896248240553  0.9835077947659414  0.9945137498637422
0.0047282952713914  0.4827940046010823  0.4942929782542009
0.5040349492352973  0.9821499237428384  0.4944941755970405
0.5058945352747628  0.4828206016032297  0.9940420309511140
0.5045613848356609  0.4811103128264023  0.4933877756337353
0.0013796816599694  0.9829379087234287  0.4953360299212051
0.0083465288988691  0.4848714537370853  0.9941122597789658""")
    actual = symmetrize_defect_structure(structure=structure,
                                         anchor_atom_idx=1,
                                         anchor_atom_coord=[0.0, 0.5, 0.5])
    expected = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
5 0 0
0 5 0
0 0 5
Mg   O
4     3
Direct
0.0 0.0 0.0
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
0.5 0.5 0.5
0.0 0.0 0.5
0.0 0.5 0.0""")
    assert actual == expected


"""
TODO
- Add anchor_atom_index

DONE
"""