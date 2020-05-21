# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip


def test_calc_results_msonable(vasp_files, tmpdir):
    vasprun = Vasprun(vasp_files / "MgO_conv_Va_O_0" / "vasprun.xml")
    outcar = Outcar(vasp_files / "MgO_conv_Va_O_0" / "OUTCAR")
    calc_results = make_calc_results_from_vasp(vasprun, outcar)
    assert_msonable(calc_results)
    assert_json_roundtrip(calc_results, tmpdir)

