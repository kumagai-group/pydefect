# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
from monty.serialization import loadfn
from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.vasp.make_gkfo_correction import make_gkfo_correction
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection
from pydefect.corrections.gkfo_correction import GkfoCorrection


def test_make_gkfo_correction(vasp_files):
    d = vasp_files / "MgO_VaO_+_to_2+_transition"
    efnv_corr: ExtendedFnvCorrection = loadfn(d / "initial_+" / "correction.json")
    unitcell: Unitcell = loadfn(d / "unitcell.json")
    initial_calc_results = loadfn(d / "initial_+" / "calc_results.json")
    final_calc_results = loadfn(d / "final_2+" / "calc_results.json")

    gkfo = make_gkfo_correction(
        efnv_correction=efnv_corr,
        additional_charge=1,
        final_calc_results=final_calc_results,
        initial_calc_results=initial_calc_results,
        diele_tensor=unitcell.dielectric_constant,
        ion_clamped_diele_tensor=unitcell.ele_dielectric_const)

    np.testing.assert_almost_equal(gkfo.pc_1st_term, 0.5578636006366229)
    np.testing.assert_almost_equal(gkfo.pc_2nd_term, 0.883976881189451)
    np.testing.assert_almost_equal(gkfo.alignment_1st_term, -0.0743015221926139)
    np.testing.assert_almost_equal(gkfo.alignment_2nd_term, -0.1299646963025651)
    np.testing.assert_almost_equal(gkfo.alignment_3rd_term, -0.0234452481649620)
    np.testing.assert_almost_equal(gkfo.correction_energy, 1.2141290151659327)

    gkfo.to_json_file("gkfo_correction.json")




