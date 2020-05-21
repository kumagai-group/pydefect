# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pymatgen import Structure
from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp


def test_make_calc_results_from_vasp_results(vasp_files):
    vasprun = Vasprun(vasp_files / "MgO_conv_Va_O_0" / "vasprun.xml")
    outcar = Outcar(vasp_files / "MgO_conv_Va_O_0" / "OUTCAR")
    results = make_calc_results_from_vasp(vasprun, outcar)

    expected_structure = \
        Structure.from_file(vasp_files / "MgO_conv_Va_O_0" / "CONTCAR")
    assert results.structure == expected_structure
    assert results.site_symmetry == "m-3m"
    assert results.energy == -34.91084360
    assert results.magnetization == 1.03e-05

