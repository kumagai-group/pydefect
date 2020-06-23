# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import numpy as np
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues
from pydefect.defaults import defaults
from pymatgen.io.vasp import Vasprun


def make_band_edge_eigenvalues(vasprun: Vasprun, vbm: float, cbm: float):
    """

    eigenvalues[Spin][k-index][band-index] = [energy, occupation]
    """
    eigval_range = defaults.eigval_range
    kpt_coords = [tuple(coord) for coord in vasprun.actual_kpoints]

    max_energy_by_spin = []
    min_energy_by_spin = []

    for e in vasprun.eigenvalues.values():
        max_energy_by_spin.append(np.amax(e[:, :, 0], axis=0))
        min_energy_by_spin.append(np.amin(e[:, :, 0], axis=0))

    max_energy_by_band = np.amax(np.vstack(max_energy_by_spin), axis=0)
    min_energy_by_band = np.amin(np.vstack(min_energy_by_spin), axis=0)

    lower_idx = np.argwhere(max_energy_by_band > vbm - eigval_range)[0][0]
    upper_idx = np.argwhere(min_energy_by_band < cbm + eigval_range)[-1][-1]

    e_and_occu = []
    for e in vasprun.eigenvalues.values():
        e_and_occu.append(e[:, lower_idx:upper_idx + 1, :].tolist())

    return BandEdgeEigenvalues(energies_and_occupations=e_and_occu,
                               kpt_coords=kpt_coords,
                               # need to convert numpy.int64 to int for mongoDB.
                               lowest_band_index=int(lower_idx))
