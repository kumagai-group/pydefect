# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.chem_pot_diag.datasets.make_atom_vasp_set import \
    make_atom_mp_relax_set


#
#
# def test(tmpdir, mocker):
#     print(tmpdir)
#     tmpdir.chdir()
#     a = PotcarSet.normal
#     make_atom_vasp_set(a, Xc.pbesol)


def test_atom_mp_relax_set(tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    make_atom_mp_relax_set()


"""
TODO
- Make H dir with potcar inside

DONE
"""
