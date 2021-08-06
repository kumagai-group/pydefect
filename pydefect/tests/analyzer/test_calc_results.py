# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from pydefect.analyzer.calc_results import CalcResults
from pymatgen.core import IStructure, Lattice

from vise.tests.helpers.assertion import assert_json_roundtrip


@pytest.fixture
def calc_results():
    return CalcResults(structure=IStructure(Lattice.cubic(1.0), ["H"], [[0.0]*3]),
                       energy=1.0,
                       magnetization=0.0,
                       potentials=[0.0],
                       electronic_conv=False,
                       ionic_conv=False)


def test_calc_results_json_roundtrip(calc_results, tmpdir):
    assert_json_roundtrip(calc_results, tmpdir)


def test_calc_results_show_convergence_warning(calc_results):
    calc_results.show_convergence_warning()


def test_calc_results_str(calc_results):
    assert calc_results.__str__() == """ -- calc results info
energy:      1.000
magnetization:   0.00
electronic convergence: False
ionic convergence: False"""
