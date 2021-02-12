# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest

from pathlib import Path

from pymatgen import Lattice, DummySpecies
from pymatgen.core.structure import Structure

from pydefect.analyzer.vesta.create_vesta_file import VestaFile, Title, Cellp, Struc, Bound, SBond, Vect, Style

parent_dir = Path(__file__).parent


@pytest.fixture
def structure():
    return Structure(Lattice.cubic(3.9), species=["Ba", DummySpecies()], coords=[[0.5]*3, [0.0]*3])


def test_title(structure):
    actual = repr(Title(structure))
    expected = '''#VESTA_FORMAT_VERSION 3.5.0
TITLE
X1 Ba1'''
    assert actual == expected


def test_cellp(structure):
    actual = repr(Cellp(structure))
    expected = '''CELLP
3.900000 3.900000 3.900000 90.000000 90.000000 90.000000'''
    assert actual == expected


def test_struc(structure):
    actual = repr(Struc(structure))
    expected = '''STRUC
1 Ba Ba1 1.0 0.500000 0.500000 0.500000
 0.0 0.0 0.0 
2 XX XX2 1.0 0.000000 0.000000 0.000000
 0.0 0.0 0.0 
 0 0 0 0 0 '''
    assert actual == expected


def test_bound():
    actual = repr(Bound((0, 2, 0, 2, 0, 2)))
    expected = '''BOUND
0.000000 2.000000 0.000000 2.000000 0.000000 2.000000
 0 0 0 0 0 '''
    assert actual == expected


def test_sbond():
    s = Structure(Lattice.cubic(3.9), species=["Ba", "Ti"], coords=[[0.5]*3, [0.0]*3])
    actual = repr(SBond(s))
    expected = '''SBOND
Ba Ti 0.0  2.81  0  1  1  0  1
Ti Ba 0.0  2.81  0  1  1  0  1
 0 0 0 0 '''
    assert actual == expected


def test_vect():
    actual = repr(Vect({1: [0., 0., 0.], 2: [0., 0., -0.1]}))
    expected = '''VECTR
1 0.000000 0.000000 0.000000
1  0 0 0 0  0 0 0 0 0 
 0 0 0 0 0 
2 0.000000 0.000000 -0.100000
2  0 0 0 0  0 0 0 0 0 
 0 0 0 0 0 

VECTT
1 0.5 1 1 1 2
2 0.5 1 1 1 2'''
    assert actual == expected


def test_style(structure):
    actual = repr(Style(structure, is_ionic=True))
    expected = '''STYLE
VECTS  9.8865
SECTS  160  1
SECTP 
   1 0 0  0.00000E+00  0.00000E+00  0.00000E+00  0.00100E+00
UCOLP 
 0  2  1.000   0   0   0
ATOMS  1  0  1'''
    assert actual == expected


def test_vesta_file(structure):
    print(VestaFile(structure))

