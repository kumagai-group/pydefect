# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_analyzer import \
    symmetrize_defect_structure, \
    DefectStructureAnalyzer
from pymatgen import Structure, IStructure, Lattice


@pytest.fixture
def structure_analyzer():
    cu2o_perfect = IStructure(Lattice.cubic(5),
                              species=["Cu"] * 4 + ["O"] * 2,
                              coords=[[0.25, 0.25, 0.25],
                                      [0.25, 0.75, 0.75],
                                      [0.75, 0.75, 0.25],
                                      [0.75, 0.25, 0.75],
                                      [0, 0, 0],
                                      [0.5, 0.5, 0.5]])
    cu2o_defect = IStructure(Lattice.cubic(5),
                             species=["Cu"] * 3 + ["O"] * 2,
                             coords=[[0.25, 0.5, 0.5],
                                     [0.76, 0.73, 0.24],
                                     [0.75, 0.25, 0.73],
                                     [0.1, -0.1, 0],
                                     [0.5, 0.5, 0.5]])

    return DefectStructureAnalyzer(defective_structure=cu2o_defect,
                                   perfect_structure=cu2o_perfect)


def test_atom_mapping_to_perfect(structure_analyzer):
    assert structure_analyzer.atom_mapping == {1: 2, 2: 3, 3: 4, 4: 5}
    assert structure_analyzer.vacancy_indices == [0, 1]
    assert structure_analyzer.inserted_indices == [0]


def test_defect_structure_analyzer_defect_center(structure_analyzer):
    expected = np.array([0.25, 0.5, 0.5])
    assert (structure_analyzer.defect_center_coord == expected).all()


def test_defect_structure_analyzer_distance(structure_analyzer):
    assert structure_analyzer.distance_from_center(4) == 5 / 4


def test_neighboring_atom_indices(cubic_supercell):
    structure: Structure = cubic_supercell.copy()
    structure.pop(32)  # [0.25, 0, 0]
    structure.pop(0)  # [0, 0, 0]
    structure.append(species="Li", coords=[0.124, 0, 0])
    structure.to(filename="POSCAR")
    structure_analyzer = DefectStructureAnalyzer(structure, cubic_supercell)
    assert structure_analyzer.neighboring_atom_indices == sorted([25, 15, 16, 23, 54, 47, 46, 56, 62])
    # must be list to be used for indexing


def test_actual_files(vasp_files):
    d = vasp_files / "KInO2_Va_O_2"
    calc_results: CalcResults = loadfn(d / "Va_O1_2/calc_results.json")
    perfect_calc_results: CalcResults = loadfn(d / "perfect_calc_results.json")
    dsa = DefectStructureAnalyzer(calc_results.structure,
                                  perfect_calc_results.structure)
    a = list(range(192))
    a.pop(96)
    expected = dict(zip(range(191), a))
    assert dsa.inserted_indices == []
    assert dsa.vacancy_indices == [96]
    assert dsa.atom_mapping == expected


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