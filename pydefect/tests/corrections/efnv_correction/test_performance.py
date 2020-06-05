# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest
from monty.serialization import loadfn

from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.vasp.make_efnv_correction import \
    make_efnv_correction
from pydefect.input_maker.defect_entry import DefectEntry


@pytest.mark.skip()
def test(vasp_files):
    d = vasp_files / "KInO2_Va_O_2"
    defect_entry: DefectEntry = loadfn(d / "Va_O1_2/defect_entry.json")
    calc_results: CalcResults = loadfn(d / "Va_O1_2/calc_results.json")
    perfect_calc_results: CalcResults = loadfn(d / "perfect_calc_results.json")
    unitcell: Unitcell = loadfn(d / "unitcell.json")

    make_efnv_correction(defect_entry.charge,
                         calc_results,
                         perfect_calc_results,
                         unitcell.dielectric_constant,
                         1)


