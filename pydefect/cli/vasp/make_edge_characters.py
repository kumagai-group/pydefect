# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from typing import Dict

import numpy as np
from pydefect.analyzer.band_edge_states import EdgeCharacters, EdgeCharacter
from pymatgen import Spin, Structure
from pymatgen.io.vasp import Procar, Vasprun, Outcar
from vise.analyzer.vasp.band_edge_properties import eigenvalues_from_vasprun


class MakeEdgeCharacters:
    def __init__(self,
                 procar: Procar,
                 vasprun: Vasprun,
                 outcar: Outcar,
                 neighboring_atom_indices):
        self.orbs = procar.data
        self.structure = vasprun.final_structure
        self.eigenvalues = eigenvalues_from_vasprun(vasprun)
        self.nelect = outcar.nelect
        self.mag = outcar.total_mag if outcar.total_mag else 0.0
        self.neighboring_atom_indices = neighboring_atom_indices

    @property
    def edge_characters(self):
        chars = []
        for spin, eigenvalues in self.eigenvalues.items():
            if spin == Spin.up:
                vbm_band_idx = int(round((self.nelect + self.mag) / 2)) - 1
            else:
                vbm_band_idx = int(round((self.nelect - self.mag) / 2)) - 1
            cbm_band_idx = vbm_band_idx + 1

            vbm = np.amax(eigenvalues[:, vbm_band_idx], axis=0)
            vbm_kpt_idx = np.argwhere(eigenvalues[:, vbm_band_idx] == vbm)[0][0]
            hob_bottom_e = np.amin(eigenvalues[:, vbm_band_idx], axis=0)

            cbm = np.amin(eigenvalues[:, cbm_band_idx], axis=0)
            cbm_kpt_idx = np.argwhere(eigenvalues[:, cbm_band_idx] == cbm)[0][0]
            lub_top_e = np.amax(eigenvalues[:, cbm_band_idx], axis=0)

            vbm_character = calc_orbital_character(self.orbs, self.structure, spin, vbm_kpt_idx, vbm_band_idx)
            cbm_character = calc_orbital_character(self.orbs, self.structure, spin, cbm_kpt_idx, cbm_band_idx)

            if self.neighboring_atom_indices:
                vbm_participation_ratio = calc_participation_ratio(self.orbs, spin, vbm_kpt_idx, vbm_band_idx, self.neighboring_atom_indices)
                cbm_participation_ratio = calc_participation_ratio(self.orbs, spin, cbm_kpt_idx, cbm_band_idx, self.neighboring_atom_indices)
            else:
                vbm_participation_ratio = None
                cbm_participation_ratio = None

            chars.append(EdgeCharacter(hob_bottom_e, lub_top_e, vbm, cbm, vbm_character, cbm_character, vbm_participation_ratio, cbm_participation_ratio))

        return EdgeCharacters(chars)


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




