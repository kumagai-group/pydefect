# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections import defaultdict
from itertools import zip_longest

from pydefect.analyzer.band_edge_states import EdgeState, EdgeCharacter
from pydefect.defaults import defaults


def make_band_edge_state(target: EdgeCharacter, ref: EdgeCharacter) -> EdgeState:
    is_hob_near_cbm = (target.hob_bottom_e
                       > ref.cbm - defaults.similar_energy_criterion)
    is_lub_near_vbm = (target.lub_top_e
                       < ref.vbm + defaults.similar_energy_criterion)

    if target.vbm is not None:
        similar_vbm_e = (abs(target.vbm - ref.vbm)
                         < defaults.similar_energy_criterion)
    else:
        similar_vbm_e = False

    if target.cbm is not None:
        similar_cbm_e = (abs(target.cbm - ref.cbm)
                         < defaults.similar_energy_criterion)
    else:
        similar_cbm_e = False

    similar_vbm_o = are_orbitals_similar(target.vbm_orbitals, ref.vbm_orbitals)
    similar_cbm_o = are_orbitals_similar(target.cbm_orbitals, ref.cbm_orbitals)

    if similar_vbm_e and similar_cbm_e and similar_vbm_o and similar_cbm_o:
        return EdgeState.no_in_gap
    elif is_hob_near_cbm:
        return EdgeState.donor_phs
    elif is_lub_near_vbm:
        return EdgeState.acceptor_phs
    else:
        return EdgeState.in_gap_state


def are_orbitals_similar(orbital_1: dict, orbital_2: dict) -> float:
    element_set = set(list(orbital_1.keys()) + list(orbital_2.keys()))
    orb_1, orb_2 = defaultdict(list, orbital_1), defaultdict(list, orbital_2)

    difference = 0
    for e in element_set:
        difference += sum([abs(i - j) for i, j
                           in zip_longest(orb_1[e], orb_2[e], fillvalue=0)])

    return difference < defaults.similar_orb_criterion
