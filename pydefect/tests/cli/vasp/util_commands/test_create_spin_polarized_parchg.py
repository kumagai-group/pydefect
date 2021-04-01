# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.cli.vasp.util_commands.create_spin_polarized_parchg import \
    make_spin_decomposed_parchg
from pymatgen import Structure, Lattice
from pymatgen.io.vasp import Chgcar
import numpy as np


def test_make_spin_decomposed_parchg(tmpdir):
    tmpdir.chdir()
    structure = Structure(Lattice.cubic(10.0), species=["H"], coords=[[0.0]*3])
    data = {"total": np.array([[[1.1, 0.0]]]), "diff": np.array([[[0.2, 0.0]]])}
    parchg = Chgcar(structure, data)
    make_spin_decomposed_parchg(parchg, 10)
    parchg_up = Chgcar.from_file("PARCHG_10_up")
    parchg_down = Chgcar.from_file("PARCHG_10_down")
    np.testing.assert_array_almost_equal(parchg_up.data["total"], np.array([[[0.65, 0.0]]]))
    np.testing.assert_array_almost_equal(parchg_down.data["total"], np.array([[[0.45, 0.0]]]))


def test_make_spin_decomposed_parchg_2(tmpdir):
    tmpdir.chdir()
    print(tmpdir)
    structure = Structure(Lattice.cubic(10.0), species=["H"], coords=[[0.0]*3])
    data = {"total": np.array([[[1.1, 0.0]]]), "diff": np.array([[[0.2, 0.0]]])}
    parchg = Chgcar(structure, data)
    make_spin_decomposed_parchg(parchg, 10, [0.5])
    parchg_up = Chgcar.from_file("PARCHG_10_up")
    parchg_down = Chgcar.from_file("PARCHG_10_down")
    np.testing.assert_array_almost_equal(parchg_up.data["total"], np.array([[[1, 0]]]))
    np.testing.assert_array_almost_equal(parchg_down.data["total"], np.array([[[1, 0]]]))


"""
- Add info to the header
- compress option
- Receive defect.vesta and PARCHG_10_up, and create defect_structure_10_up.vesta
-
"""