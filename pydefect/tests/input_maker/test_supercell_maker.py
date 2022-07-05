# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
import pytest
from pydefect.input_maker.supercell_maker import SupercellMaker
from pydefect.util.error_classes import NotPrimitiveError
from pymatgen.core import Lattice, Structure
from vise.tests.helpers.assertion import assert_structure_almost_same


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
        SupercellMaker(primitive_structure=bcc)


def test_create_supercell_matrix(a_centered_orthorhombic):
    matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    cs = SupercellMaker(a_centered_orthorhombic, matrix_to_conv_cell=matrix)
    actual = cs.supercell.matrix
    np.testing.assert_array_equal(actual, matrix)

    expected = [[2, 0, 0], [0, 1, -1], [0, 1, 1]]
    np.testing.assert_array_equal(cs.transformation_matrix, expected)


def test_create_supercell_generate_supercell_info(simple_cubic):
    cs = SupercellMaker(simple_cubic, matrix_to_conv_cell=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
    info = cs.supercell_info
    assert info.transformation_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    assert info.sites["H1"].wyckoff_letter == "a"
    assert info.sites["H1"].equivalent_atoms == list(range(8))
    assert info.unitcell_structure == simple_cubic


def test_spglib_cyclic_behavior():
    structure = Structure.from_str("""Ca4 Sc2 Sb2 O12
1.0
5.469740 0.000000 0.000000
0.000000 5.632621 0.000000
-5.471566 0.000000 7.830251
O
12
direct
   0.242539  0.802107  0.449144
   0.242539  0.697893  0.949144
   0.757461  0.197893  0.550856
   0.757461  0.302107  0.050856
   0.859254  0.791326  0.05328
   0.859254  0.708674  0.55328
   0.140746  0.208674  0.94672
   0.140746  0.291326  0.44672
   0.656928  0.965329  0.256023
   0.656928  0.534671  0.756023
   0.343072  0.034671  0.743977
   0.343072  0.465329  0.243977
""", fmt="poscar")

    sm = SupercellMaker(structure)
    assert sm.conv_structure == structure


def test_transformation_matrix():
    structure = Structure.from_str("""
   1.0
    -2.5007160000000002    2.5007160000000002    3.6755999999999993
     2.5007160000000002   -2.5007159999999988    3.6755999999999993
     2.5007160000000002    2.5007159999999988   -3.6755999999999993
Si
   2    
Direct
  0.7500000000000001  0.2500000000000001  0.5000000000000000
  0.0000000000000000  0.0000000000000000  0.0000000000000000""", fmt="poscar")

    conv_structure = Structure.from_str("""
  1.0
5.001432 0.000000 0.000000
0.000000 5.001432 0.000000
0.000000 0.000000 7.351200
Si 
4 
direct
0.000000 0.500000 0.250000 Si
0.500000 1.000000 0.750000 Si
0.000000 0.000000 0.000000 Si
0.500000 0.500000 0.500000 Si""", fmt="poscar")

    matrix_to_conv_cell = [[3, 0, 0], [0, 3, 0], [0, 0, 2]]
    sm = SupercellMaker(structure,
                        matrix_to_conv_cell=matrix_to_conv_cell,
                        symprec=0.01)
    actual = structure * sm.transformation_matrix
    expected = conv_structure * matrix_to_conv_cell
    assert_structure_almost_same(actual, expected)


def test_rhombohedral_supercell_info():
    """There is a bug in spglib to find the inequivalent sites in supercell"""
    structure = Structure.from_str("""    Ba4 Nb2 W1 O12
1.00000000000000
2.8901645462281422    1.6686372788168120    9.3650384663595574
-2.8901645462281422    1.6686372788168120    9.3650384663595574
0.0000000000000000   -3.3372745575340472    9.3650384663596231
Ba   Nb   W    O
4     2     1    12
Direct
0.2900069010571258  0.2900069010571258  0.2900069010571258
0.7114602455184809  0.7114602455184809  0.7114602455184809
0.1345447049748643  0.1345447049748643  0.1345447049748643
0.8646114100059918  0.8646114100059918  0.8646114100059918
0.0006723681882974  0.0006723681882974  0.0006723681882974
0.4240665571167241  0.4240665571167241  0.4240665571167241
0.5787478829344508  0.5787478829344508  0.5787478829344508
0.7884796629179860  0.2974972818896759  0.7884796629179860
0.2974972818896759  0.7884796629179860  0.7884796629179860
0.7884796629179860  0.7884796629179860  0.2974972818896759
0.2034890710071906  0.7099376271089710  0.2034890710071906
0.7099376271089710  0.2034890710071906  0.2034890710071906
0.2034890710071906  0.2034890710071906  0.7099376271089710
0.6204119903421983  0.1229791605843715  0.6204119903421983
0.1229791605843715  0.6204119903421983  0.6204119903421983
0.6204119903421983  0.6204119903421983  0.1229791605843715
0.3875918889522652  0.8655326340817453  0.3875918889522652
0.8655326340817453  0.3875918889522652  0.3875918889522652
0.3875918889522652  0.3875918889522652  0.8655326340817453""", fmt="poscar")

    sm = SupercellMaker(structure, symprec=0.01, min_num_atoms=60, max_num_atoms=500)
    assert "O5" not in sm.supercell_info.sites.keys()