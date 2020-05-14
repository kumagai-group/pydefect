# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.util.centering import Centering
from vise.util.structure_symmetrizer import StructureSymmetrizer


def test_names():
    assert Centering.names_string() == "P, A, C, R, I, F"


def test_a_centering(a_centered_orthorhombic):
    s = StructureSymmetrizer(a_centered_orthorhombic)
    primitive = s.primitive
    to_conventional = primitive * Centering.A.primitive_to_conv
    assert to_conventional.lattice == s.conventional.lattice
    assert Centering.A.conv_multiplicity == 2


def test_c_centering(c_centered_monoclinic):
    s = StructureSymmetrizer(c_centered_monoclinic)
    primitive = s.primitive
    to_conventional = primitive * Centering.C.primitive_to_conv
    assert to_conventional == s.conventional
    assert Centering.C.conv_multiplicity == 2


def test_r_centering(rhombohedral):
    s = StructureSymmetrizer(rhombohedral)
    rhombohedral = s.primitive
    to_conventional = rhombohedral * Centering.R.primitive_to_conv
    assert to_conventional == s.conventional
    assert Centering.R.conv_multiplicity == 3


def test_f_centering(bcc):
    s = StructureSymmetrizer(bcc)
    primitive = s.primitive
    to_conventional = primitive * Centering.I.primitive_to_conv
    assert to_conventional == s.conventional
    assert Centering.I.conv_multiplicity == 2


"""
TODO
- 

DONE
"""