# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState, EdgeInfo, \
    OrbitalInfo
from pydefect.cli.vasp.make_band_edge_orbital_infos import \
    calc_orbital_character
from vise.defaults import defaults as v_defaults
from pymatgen import Spin
from pymatgen.io.vasp import Procar, Vasprun, Outcar
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties


def make_perfect_band_edge_state_from_vasp(
        procar: Procar, vasprun: Vasprun, outcar: Outcar
) -> PerfectBandEdgeState:

    band_edge_prop = VaspBandEdgeProperties(vasprun, outcar,
                                            v_defaults.integer_criterion)
    orbs, s = procar.data, vasprun.final_structure
    vbm_info = get_edge_info(band_edge_prop.vbm_info, orbs, s, vasprun)
    cbm_info = get_edge_info(band_edge_prop.cbm_info, orbs, s, vasprun)
    return PerfectBandEdgeState(vbm_info, cbm_info)


def get_edge_info(ei, orbs, s, vasprun):
    orbitals = calc_orbital_character(
        orbs, s, Spin.up, ei.kpoint_index, ei.band_index)
    e, occ = vasprun.eigenvalues[Spin.up][ei.kpoint_index, ei.band_index, :]
    orb_info = OrbitalInfo(energy=e, occupation=occ, orbitals=orbitals)
    return EdgeInfo(ei.band_index, tuple(ei.kpoint_coords), orb_info)


