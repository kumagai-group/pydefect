# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.analyzer.calc_results import CalcResults
from pymatgen.io.vasp import Vasprun, Outcar


def make_calc_results_from_vasp(vasprun: Vasprun,
                                outcar: Outcar) -> CalcResults:
    return CalcResults(structure=vasprun.final_structure,
                       energy=outcar.final_energy,
                       magnetization=outcar.total_mag or 0.0,
                       potentials=[-p for p in outcar.electrostatic_potential],
                       electronic_conv=vasprun.converged_electronic,
                       ionic_conv=vasprun.converged_ionic)
