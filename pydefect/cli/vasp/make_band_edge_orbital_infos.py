# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import Dict

import numpy as np
from pydefect.analyzer.band_edge_states import OrbitalInfo, \
    BandEdgeOrbitalInfos
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.defaults import defaults
from pymatgen import Spin, Structure
from pymatgen.io.vasp import Procar, Vasprun


def make_band_edge_orbital_infos(procar: Procar,
                                 vasprun: Vasprun,
                                 vbm: float, cbm: float,
                                 str_info: DefectStructureInfo):
    eigval_range = defaults.eigval_range
    kpt_coords = [tuple(coord) for coord in vasprun.actual_kpoints]
    max_energy_by_spin, min_energy_by_spin = [], []
    neighbors = str_info.neighbor_atom_indices

    for e in vasprun.eigenvalues.values():
        max_energy_by_spin.append(np.amax(e[:, :, 0], axis=0))
        min_energy_by_spin.append(np.amin(e[:, :, 0], axis=0))

    max_energy_by_band = np.amax(np.vstack(max_energy_by_spin), axis=0)
    min_energy_by_band = np.amin(np.vstack(min_energy_by_spin), axis=0)

    lower_idx = np.argwhere(max_energy_by_band > vbm - eigval_range)[0][0]
    upper_idx = np.argwhere(min_energy_by_band < cbm + eigval_range)[-1][-1]

    orbs, s = procar.data, vasprun.final_structure
    orb_infos = []
    for spin, eigvals in vasprun.eigenvalues.items():
        orb_infos.append([])
        for k_idx in range(len(kpt_coords)):
            orb_infos[-1].append([])
            for b_idx in range(lower_idx, upper_idx + 1):
                e, occ = eigvals[k_idx, b_idx, :]
                orbitals = calc_orbital_character(orbs, s, spin, k_idx, b_idx)
                p_ratio = calc_participation_ratio(orbs, spin, k_idx, b_idx, neighbors)
                orb_infos[-1][-1].append(OrbitalInfo(e, orbitals, occ, p_ratio))

    return BandEdgeOrbitalInfos(
        orbital_infos=orb_infos,
        kpt_coords=kpt_coords,
        kpt_weights=vasprun.actual_kpoints_weights,
        # need to convert numpy.int64 to int for mongoDB.
        lowest_band_index=int(lower_idx),
        fermi_level=vasprun.efermi)


def calc_participation_ratio(orbitals: Dict[Spin, np.ndarray],
                             spin: Spin,
                             kpt_index: int,
                             band_index: int,
                             atom_indices: list) -> float:
    """ Returns sum of participation ratios at atom_indices sites

    The PROCAR data of the form below. It should VASP uses 1-based indexing,
    but all indices are converted to 0-based here.::
        { spin: np.ndarray accessed with (k-point index, band index,
                                          ion index, orbital index) }

    Note that the k-point weight is not considered, so all the k-points are
    treated equally.

    Return (float):
        float of the participation ratio.
    """
    sum_per_atom = np.sum(orbitals[spin][kpt_index, band_index, :, :], axis=1)
    return np.sum(sum_per_atom[atom_indices]) / np.sum(sum_per_atom)


def calc_orbital_character(orbitals,
                           structure: Structure,
                           spin: Spin,
                           kpt_index: int,
                           band_index: int):
    """ Consider the two pattern of orbitals.

    LORBIT 10 -> consider only "s", "p", "d" orbitals
    LORBIT >=11 -> consider "s", "px", "py", "pz",.. orbitals

    """

    def projection_sum(atom_indices: tuple, first: int, last: int):
        end = last + 1
        procar_sum = np.sum(orbitals[spin]
                            [kpt_index, band_index, atom_indices, first:end])
        return float(procar_sum)

    orbital_components = {}
    azimuthal = len(orbitals[Spin.up][0, 0, 0]) > 5

    for element in structure.symbol_set:
        # get list of index
        indices = structure.indices_from_symbol(element)
        if azimuthal:
            orbital_components[element] = \
                [round(projection_sum(indices, 0, 0), 3),
                 round(projection_sum(indices, 1, 3), 3),
                 round(projection_sum(indices, 4, 8), 3)]
            try:
                orbital_components[element].append(round(projection_sum(indices, 9, 16), 3))
            except KeyError:
                pass
        else:
            orbital_components[element] = \
                [round(projection_sum(indices, 0, 0), 3),
                 round(projection_sum(indices, 1, 1), 3),
                 round(projection_sum(indices, 2, 2), 3)]
    return orbital_components