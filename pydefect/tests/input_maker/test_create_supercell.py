# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest
from pymatgen import Lattice

from pydefect.input_maker.create_supercell import CreateSupercell
from pydefect.util.error_classes import NotPrimitiveError


def test_create_supercell_tetragonal(tetra_close_to_cubic):
    cs = CreateSupercell(tetra_close_to_cubic)
    actual = cs.supercell.matrix
    expected =[[ 5,  5,  0],
               [-5,  5,  0],
               [ 0,  0,  5]]
    np.testing.assert_array_equal(actual, expected)


def test_create_supercell(a_centered_orthorhombic):
    cs = CreateSupercell(a_centered_orthorhombic)
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
        CreateSupercell(input_structure=bcc)


def test_create_supercell_matrix(a_centered_orthorhombic):
    matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    cs = CreateSupercell(a_centered_orthorhombic, matrix=matrix)
    actual = cs.supercell.matrix
    np.testing.assert_array_equal(actual, matrix)

    expected = [[2, 0, 0], [0, 1, -1], [0, 1, 1]]
    np.testing.assert_array_equal(cs.transform_matrix, expected)


def test_create_supercell_generate_supercell_info(simple_cubic):
    cs = CreateSupercell(simple_cubic, matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    info = cs.supercell_info
    assert info.transform_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    assert info.sites["H1"].wyckoff_letter == "a"
    assert info.sites["H1"].equivalent_atoms == list(range(8))

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