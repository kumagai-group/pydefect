# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Set

import numpy as np
from monty.json import MSONable
from pydefect.defaults import defaults
from pydefect.util.coords import pretty_coords
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.typing import Coords, GenCoords

printed_orbital_weight_threshold = 0.1


@dataclass
class BandEdgeEigenvalues(MSONable, ToJsonFileMixIn):
    # [spin, k-idx, band-idx] = energy, occupation
    energies_and_occupations: List[List[List[List[float]]]]
    kpt_coords: List[List[float]]
    lowest_band_index: int


def pretty_orbital(orbitals: Dict[str, List[float]]):
    """
    :param orbitals: An example is {"Mn": [0.5, 0.4, 0.01]} (no f-orbital)
    :return: "Mn-s: 0.50, Mn-p: 0.40"
    """
    orbital_infos = []
    for elem, orbs in orbitals.items():
        for orb_name, weight in zip(["s", "p", "d", "f"], orbs):
            if weight > printed_orbital_weight_threshold:
                orbital_infos.append(f"{elem}-{orb_name}: {weight:.2f}")
    return ", ".join(orbital_infos)


@dataclass
class OrbitalInfo(MSONable):
    """Note that this is the code- and version-dependent quantities. """
    energy: float  # max eigenvalue
    # {"Mn": [0.01, ..], "O": [0.03, 0.5]},
    # where lists contain s, p, d, (f) orbital components.
    orbitals: Dict[str, List[float]]
    occupation: float
    participation_ratio: float = None


@dataclass
class BandEdgeOrbitalInfos(MSONable, ToJsonFileMixIn):
    orbital_infos: List[List[List["OrbitalInfo"]]]  # [spin, k-idx, band-idx]
    kpt_coords: List[Coords]
    kpt_weights: List[float]
    lowest_band_index: int
    fermi_level: float

    def kpt_idx(self, kpt_coord):
        for i, orig_kpt in enumerate(self.kpt_coords):
            if sum(abs(k - l) for k, l in zip(orig_kpt, kpt_coord)) < 1e-5:
                return i
        raise ValueError(f"{kpt_coord} not in kpt_coords {self.kpt_coords}.")

    @property
    def energies_and_occupations(self) -> List[List[List[List[float]]]]:
        result = np.zeros(np.shape(self.orbital_infos) + (2,))
        for i, x in enumerate(self.orbital_infos):
            for j, y in enumerate(x):
                for k, z in enumerate(y):
                    result[i][j][k] = [z.energy, z.occupation]
        return result.tolist()

    def __str__(self):
        return "\n".join([" -- band-edge orbitals info",
                          "K-points info",
                          self._kpt_block,
                          "",
                          "Band info near band edges",
                          self._band_block])

    @property
    def _band_block(self):
        band_block = [["Index", "Kpoint index", "Energy", "Occupation",
                       "P-ratio", "Orbital"]]
        for orbital_info in self.orbital_infos:
            max_idx, min_idx, t_orb_info = self._band_idx_range(orbital_info)
            for band_idx in range(min_idx, max_idx):
                actual_band_idx = band_idx + self.lowest_band_index + 1
                for kpt_idx, orb_info in enumerate(t_orb_info[band_idx], 1):

                    energy = f"{orb_info.energy :5.2f}"
                    occupation = f"{orb_info.occupation:4.1f}"
                    if orb_info.participation_ratio:
                        p_ratio = f"{orb_info.participation_ratio:4.1f}"
                    else:
                        p_ratio = "N.A."
                    orbs = pretty_orbital(orb_info.orbitals)
                    band_block.append([actual_band_idx, kpt_idx, energy,
                                       occupation, p_ratio, orbs])
                band_block.append(["--"])
            band_block.append("")
        return tabulate(band_block, tablefmt="plain")

    @property
    def _kpt_block(self):
        kpt_block = [["Index", "Coords", "Weight"]]
        for index, (kpt_coord, kpt_weight) in enumerate(
                zip(self.kpt_coords, self.kpt_weights), 1):
            coord = pretty_coords(kpt_coord)
            weight = f"{kpt_weight:4.3f}"
            kpt_block.append([index, coord, weight])
        return tabulate(kpt_block, tablefmt="plain")

    @staticmethod
    def _band_idx_range(orbital_info: List[List[OrbitalInfo]]
                        ) -> Tuple[int, int, List[List[OrbitalInfo]]]:
        # swap [kpt_idx][band_idx] -> [bane_idx][kpt_idx]
        t_orbital_info = np.array(orbital_info).T.tolist()
        middle_idx = int(len(t_orbital_info) / 2)
        for band_idx, (upper, lower) in enumerate(zip(t_orbital_info[1:],
                                                      t_orbital_info[:-1])):
            # determine the band_idx where the occupation changes largely.
            if lower[0].occupation - upper[0].occupation > 0.1:
                middle_idx = band_idx + 1
                break
        max_idx = min(middle_idx + 3, len(t_orbital_info))
        min_idx = max(middle_idx - 3, 0)
        return max_idx, min_idx, t_orbital_info


@dataclass
class LocalizedOrbital(MSONable):
    band_idx: int
    ave_energy: float
    occupation: float
    orbitals: Dict[str, List[float]]
    participation_ratio: Optional[float] = None
    radius: Optional[float] = None
    center: Optional[GenCoords] = None


@dataclass
class EdgeInfo(MSONable):
    """Information for a particular band edge, namely VBM or CBM."""
    band_idx: int
    kpt_coord: Coords
    orbital_info: "OrbitalInfo"

    @property
    def orbitals(self):
        return self.orbital_info.orbitals

    @property
    def energy(self):
        return self.orbital_info.energy

    @property
    def occupation(self):
        return self.orbital_info.occupation

    @property
    def p_ratio(self):
        return self.orbital_info.participation_ratio


@dataclass
class PerfectBandEdgeState(MSONable, ToJsonFileMixIn):
    vbm_info: EdgeInfo
    cbm_info: EdgeInfo

    def __str__(self):
        def show_edge_info(edge_info: EdgeInfo):
            return [edge_info.band_idx,
                    edge_info.energy,
                    f"{edge_info.occupation:5.2f}",
                    pretty_orbital(edge_info.orbital_info.orbitals),
                    pretty_coords(edge_info.kpt_coord)]

        return tabulate([
            ["", "Index", "Energy", "Occupation", "Orbitals", "K-point coords"],
            ["VBM"] + show_edge_info(self.vbm_info),
            ["CBM"] + show_edge_info(self.cbm_info)], tablefmt="plain")


@dataclass
class BandEdgeState(MSONable):
    """Information related to the band edge in a defect supercell.

    E.g.,
    # 124   --   <- CBM
    # 123   --   <- localized orbital
    # 122   --   <- localized orbital
    # 121   --   <- VBM
    """

    vbm_info: EdgeInfo
    cbm_info: EdgeInfo
    vbm_orbital_diff: float  # Difference from those of perfect supercell.
    cbm_orbital_diff: float
    localized_orbitals: List[LocalizedOrbital]
    vbm_hole_occupation: float = None
    cbm_electron_occupation: float = None

    @property
    def is_shallow(self):
        return self.has_donor_phs or self.has_acceptor_phs

    @property
    def has_donor_phs(self):
        return self.cbm_electron_occupation > defaults.state_occupied_threshold

    @property
    def has_acceptor_phs(self):
        return self.vbm_hole_occupation > defaults.state_occupied_threshold

    @property
    def has_unoccupied_localized_state(self):
        return any([lo.occupation < defaults.state_unoccupied_threshold
                    for lo in self.localized_orbitals])

    @property
    def has_occupied_localized_state(self):
        return any([lo.occupation > defaults.state_occupied_threshold
                    for lo in self.localized_orbitals])

    def __str__(self):
        return "\n".join([self._edge_info,
                          "---", "Localized Orbital(s)",
                          self._orbital_info])

    @property
    def _orbital_info(self):
        inner_table = [["Index", "Energy", "P-ratio", "Occupation", "Orbitals"]]
        w_radius = self.localized_orbitals and self.localized_orbitals[0].radius
        if w_radius:
            inner_table[0].extend(["Radius", "Center"])

        for lo in self.localized_orbitals:
            participation_ratio = f"{lo.participation_ratio:5.2f}" \
                if lo.participation_ratio else "None"
            inner = [lo.band_idx + 1,
                     f"{lo.ave_energy:7.3f}",
                     participation_ratio,
                     f"{lo.occupation:5.2f}",
                     pretty_orbital(lo.orbitals)]
            if w_radius:
                inner.extend([f"{lo.radius:5.2f}", pretty_coords(lo.center)])
            inner_table.append(inner)

        return tabulate(inner_table, tablefmt="plain")

    @property
    def _edge_info(self):
        inner_table = [["", "Index", "Energy", "P-ratio", "Occupation",
                        "OrbDiff", "Orbitals", "K-point coords"],
                       ["VBM"] + self._show_edge_info(
                           self.vbm_info, self.vbm_orbital_diff),
                       ["CBM"] + self._show_edge_info(
                           self.cbm_info, self.cbm_orbital_diff)]
        table = tabulate(inner_table, tablefmt="plain")
        vbm_phs = f"vbm has acceptor phs: {self.has_acceptor_phs} " \
                  f"({self.vbm_hole_occupation:5.3f} vs. {defaults.state_occupied_threshold})"
        cbm_phs = f"cbm has donor phs: {self.has_donor_phs} " \
                  f"({self.cbm_electron_occupation:5.3f} vs. {defaults.state_occupied_threshold})"
        return "\n".join([table, vbm_phs, cbm_phs])

    @staticmethod
    def _show_edge_info(edge_info: EdgeInfo, orb_diff: float):
        p_ratio = f"{edge_info.p_ratio:5.2f}" if edge_info.p_ratio else "None"
        return [edge_info.band_idx + 1,
                f"{edge_info.energy:7.3f}",
                p_ratio,
                f"{edge_info.occupation:5.2f}",
                f"{orb_diff:5.2f}",
                pretty_orbital(edge_info.orbital_info.orbitals),
                pretty_coords(edge_info.kpt_coord)]


@dataclass
class BandEdgeStates(MSONable, ToJsonFileMixIn):
    states: List[BandEdgeState]  # by spin.

    @property
    def is_shallow(self):
        return any([i.is_shallow for i in self.states])

    @property
    def has_donor_phs(self):
        return any([i.has_donor_phs for i in self.states])

    @property
    def has_acceptor_phs(self):
        return any([i.has_acceptor_phs for i in self.states])

    @property
    def has_unoccupied_localized_state(self):
        return any([i.has_unoccupied_localized_state for i in self.states])

    @property
    def has_occupied_localized_state(self):
        return any([i.has_occupied_localized_state for i in self.states])

    @property
    def band_indices_from_vbm_to_cbm(self) -> List[int]:
        indices_set = set()
        for state in self.states:
            indices_set.add(state.vbm_info.band_idx)
            for lo in state.localized_orbitals:
                indices_set.add(lo.band_idx)
            indices_set.add(state.cbm_info.band_idx)
        return sorted([i for i in indices_set])

    def __str__(self):
        lines = [" -- band-edge states info"]
        for spin, state in zip(["up", "down"], self.states):
            lines.append(f"Spin-{spin}")
            lines.append(state.__str__())
            lines.append("")

        return "\n".join(lines)
