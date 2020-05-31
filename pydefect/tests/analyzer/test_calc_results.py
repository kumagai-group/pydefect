# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.make_calc_results import make_calc_results_from_vasp
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip


@pytest.fixture
def calc_results(vasp_files):
    vasprun = Vasprun(vasp_files / "MgO_conv_Va_O_0" / "vasprun.xml")
    outcar = Outcar(vasp_files / "MgO_conv_Va_O_0" / "OUTCAR")
    return make_calc_results_from_vasp(vasprun, outcar)


def test_calc_results_msonable(calc_results):
    assert_msonable(calc_results)


def test_calc_results_json_roundtrip(calc_results, tmpdir):
    assert_json_roundtrip(calc_results, tmpdir)


def test_calc_results_vbm_cbm(calc_results):
    assert calc_results.vbm == 2.4163
    assert calc_results.cbm == 10.2574

