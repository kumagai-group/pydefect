# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections import defaultdict
from itertools import zip_longest
from typing import List

import numpy as np
from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos, \
    PerfectBandEdgeState, BandEdgeStates, EdgeInfo, OrbitalInfo, \
    LocalizedOrbital, BandEdgeState
from pydefect.defaults import defaults


def get_similar_orb_idx(orbs: List[OrbitalInfo], edge_info: EdgeInfo):
    for i, orb in enumerate(orbs):
        if abs(orb.energy - edge_info.energy) < defaults.similar_energy_criterion:
            orb_diff = orbital_diff(orb.orbitals, edge_info.orbitals)
            if orb_diff < defaults.similar_orb_criterion:
                return i, orb_diff
    raise ValueError(f"Similar orbital to the edge cannot be found. "
                     f"Try to increase the similar criterion.")


def get_localized_orbs(orb_info_by_spin: List[List[OrbitalInfo]],
                       loc_band_index_range: List[int],
                       lowest_band_idx: int,
                       weights: List[float]):
    start = loc_band_index_range[0] - lowest_band_idx
    end = loc_band_index_range[1] + 1 - lowest_band_idx
    orbs_by_band_by_kpt = np.array(orb_info_by_spin)[:, start: end].T

    result = []
    for i, orbs_by_kpt in enumerate(orbs_by_band_by_kpt, loc_band_index_range[0]):
        def get_ave(key):
            return sum([getattr(orb, key) * w for orb, w in zip(orbs_by_kpt, weights)])

        orbs = defaultdict(list)
        for orb, w in zip(orbs_by_kpt, weights):
            for elem, contrb in orb.orbitals.items():
                for y, c in enumerate(contrb):
                    try:
                        orbs[elem][y] += c
                    except IndexError:
                        orbs[elem].append(0.0)
                        orbs[elem][y] += c

        result.append(LocalizedOrbital(
            band_idx=i,
            ave_energy=get_ave("energy"),
            occupation=get_ave("occupation"),
            orbitals=dict(orbs),
            participation_ratio=get_ave("participation_ratio")))
    return result


def make_band_edge_states(orbital_infos: BandEdgeOrbitalInfos,
                          p_edge_state: PerfectBandEdgeState) -> BandEdgeStates:
    vbm_k_idx = orbital_infos.kpt_coords.index(p_edge_state.vbm_info.kpt_coord)
    cbm_k_idx = orbital_infos.kpt_coords.index(p_edge_state.cbm_info.kpt_coord)

    states = []
    for orb_info_by_spin in orbital_infos.orbital_infos:
        reversed_orbs = orb_info_by_spin[vbm_k_idx][::-1]  # by band
        i, v_diff = get_similar_orb_idx(reversed_orbs, p_edge_state.vbm_info)
        idx = orbital_infos.lowest_band_index + (len(reversed_orbs) - i - 1)
        vbm_info = EdgeInfo(band_idx=idx,
                            kpt_coord=p_edge_state.vbm_info.kpt_coord,
                            orbital_info=reversed_orbs[i])

        orbs = orb_info_by_spin[cbm_k_idx]
        i, c_diff = get_similar_orb_idx(orbs, p_edge_state.cbm_info)
        idx = orbital_infos.lowest_band_index + i
        cbm_info = EdgeInfo(band_idx=idx,
                            kpt_coord=p_edge_state.cbm_info.kpt_coord,
                            orbital_info=orbs[i])
        loc_idx_range = [vbm_info.band_idx + 1, cbm_info.band_idx - 1]
        localized_orbs = get_localized_orbs(orb_info_by_spin,
                                            loc_idx_range,
                                            orbital_infos.lowest_band_index,
                                            orbital_infos.kpt_weights)

        states.append(BandEdgeState(vbm_info=vbm_info,
                                    cbm_info=cbm_info,
                                    vbm_orbital_diff=v_diff,
                                    cbm_orbital_diff=c_diff,
                                    localized_orbitals=localized_orbs))

    return BandEdgeStates(states=states)


def orbital_diff(orbital_1: dict, orbital_2: dict) -> float:
    element_set = set(list(orbital_1.keys()) + list(orbital_2.keys()))
    orb_1, orb_2 = defaultdict(list, orbital_1), defaultdict(list, orbital_2)
    result = 0
    for e in element_set:
        result += sum([abs(i - j) for i, j
                           in zip_longest(orb_1[e], orb_2[e], fillvalue=0)])

    return result


