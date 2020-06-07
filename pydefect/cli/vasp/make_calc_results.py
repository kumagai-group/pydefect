# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.analyzer.calc_results import CalcResults
from pydefect.defaults import defaults
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties
from vise.util.structure_symmetrizer import StructureSymmetrizer


def make_calc_results_from_vasp(vasprun: Vasprun,
                                outcar: Outcar) -> CalcResults:
    structure = vasprun.final_structure
    magnetization = outcar.total_mag or 0.0

    symmetrizer = StructureSymmetrizer(
        structure,
        symprec=defaults.symmetry_length_tolerance,
        angle_tolerance=defaults.symmetry_angle_tolerance,
        time_reversal=abs(magnetization) > defaults.integer_criterion)
    band_edge_prop = VaspBandEdgeProperties(
        vasprun, outcar, defaults.integer_criterion)

    return CalcResults(structure=structure,
                       site_symmetry=symmetrizer.point_group,
                       energy=outcar.final_energy,
                       magnetization=magnetization,
                       kpoint_coords=vasprun.actual_kpoints,
                       kpoint_weights=vasprun.actual_kpoints_weights,
                       potentials=[-p for p in outcar.electrostatic_potential],
                       vbm_info=band_edge_prop.vbm_info,
                       cbm_info=band_edge_prop.cbm_info,
                       fermi_level=vasprun.efermi,
                       electronic_conv=vasprun.converged_electronic,
                       ionic_conv=vasprun.converged_ionic)
