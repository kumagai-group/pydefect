# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest
from pydefect.input_maker.supercell_maker import SupercellMaker
from pydefect.util.error_classes import NotPrimitiveError
from pymatgen import Lattice, Structure


def test_create_supercell_tetragonal(tetra_close_to_cubic):
    cs = SupercellMaker(tetra_close_to_cubic)
    actual = cs.supercell.matrix
    expected =[[ 5,  5,  0],
               [-5,  5,  0],
               [ 0,  0,  5]]
    np.testing.assert_array_equal(actual, expected)


def test_create_supercell(a_centered_orthorhombic):
    cs = SupercellMaker(a_centered_orthorhombic)
    actual_lattice = cs.conv_structure.lattice
    expected = Lattice.orthorhombic(1, 4, 6)
    assert actual_lattice == expected

    actual = cs.supercell.matrix
    expected =[[6,  0,  0],
               [0,  2,  0],
               [0,  0,  1]]
    np.testing.assert_array_equal(actual, expected)


def test_create_supercell_raise_not_primitive_error(bcc):
    with pytest.raises(NotPrimitiveError):
        SupercellMaker(input_structure=bcc)


def test_create_supercell_matrix(a_centered_orthorhombic):
    matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    cs = SupercellMaker(a_centered_orthorhombic, matrix=matrix)
    actual = cs.supercell.matrix
    np.testing.assert_array_equal(actual, matrix)

    expected = [[2, 0, 0], [0, 1, -1], [0, 1, 1]]
    np.testing.assert_array_equal(cs.transformation_matrix, expected)


def test_create_supercell_generate_supercell_info(simple_cubic):
    cs = SupercellMaker(simple_cubic, matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    info = cs.supercell_info
    assert info.transformation_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    assert info.sites["H1"].wyckoff_letter == "a"
    assert info.sites["H1"].equivalent_atoms == list(range(8))


def test_spglib_cyclic_behavior():
    structure = Structure.from_str("""Ca4 Sc2 Sb2 O12
1.0
5.469740 0.000000 0.000000
0.000000 5.632621 0.000000
-5.467914 0.000000 7.830251
O
12
direct
0.344251 0.197893 0.550856 O
0.344251 0.302107 0.050856 O
0.655749 0.802107 0.449144 O
0.655749 0.697893 0.949144 O
0.752694 0.208674 0.946720 O
0.752694 0.291326 0.446720 O
0.247306 0.791326 0.053280 O
0.247306 0.708674 0.553280 O
0.144882 0.034671 0.743977 O
0.144882 0.465329 0.243977 O
0.855118 0.965329 0.256023 O
0.855118 0.534671 0.756023 O
""", fmt="POSCAR")

    sm = SupercellMaker(structure)
    assert sm.conv_structure == structure


"""
TODO


DONE
- Create SupercellInfo at CreateSupercell
- Expand for rhombohedral cell.
- Expand for tetragonal cell.
- Input of structure and trans_mat returns supercell
- Return isotropy
- Return average angle
- Recommend the least isotropic supercell of simple cubic within a given number of atom range.
  When the isotropy is the same, smallest supercell is returned.
- Check max_num_atoms
- Raise NoSupercellError 

REJECT
+ Allow to write the supercell down a file.
    -> Structure already has to, which should be used.
"""