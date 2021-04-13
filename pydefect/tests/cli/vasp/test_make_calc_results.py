# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp
from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun, Outcar


def test_make_calc_results_from_vasp_results(vasp_files):
    vasprun = Vasprun(vasp_files / "MgO_conv_Va_O_0" / "vasprun.xml")
    outcar = Outcar(vasp_files / "MgO_conv_Va_O_0" / "OUTCAR")
    results = make_calc_results_from_vasp(vasprun, outcar)

    expected_structure = \
        Structure.from_file(vasp_files / "MgO_conv_Va_O_0" / "CONTCAR")
    assert results.structure == expected_structure
    assert results.energy == -34.91084360
    assert results.magnetization == 1.03e-05
    assert results.potentials ==\
           [35.9483, 36.066, 35.948, 35.9478, 69.799, 69.7994, 69.7995]
