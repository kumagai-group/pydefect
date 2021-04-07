# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Tuple, Optional

import numpy as np
from monty.json import MSONable
from pydefect.util.coords import pretty_coords
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.typing import Coords


@dataclass
class BandEdgeEigenvalues(MSONable, ToJsonFileMixIn):
    # [spin, k-idx, band-idx] = energy, occupation
    energies_and_occupations: List[List[List[List[float]]]]
    kpt_coords: List[Tuple[float, float, float]]
    lowest_band_index: int


def pretty_orbital(orbitals):
    x = []
    for elem, orbs in orbitals.items():
        for orb_name, weight in zip(["s", "p", "d", "f"], orbs):
            if weight > 0.1:
                x.append(f"{elem}-{orb_name}: {weight:.2f}")
    return ", ".join(x)


@dataclass
class OrbitalInfo(MSONable):
    """Code and its version dependent quantities. """
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

    @property
    def energies_and_occupations(self):
        result = deepcopy(self.orbital_infos)
        for i, x in enumerate(self.orbital_infos):
            for j, y in enumerate(x):
                for k, z in enumerate(y):
                    result[i][j][k] = [z.energy, z.occupation]
        return result

    def __str__(self):
        lines = [" -- band-edge orbitals info", "K-points info"]
        kpt_block = [["Index", "Coords", "Weight"]]
        for i, (c, w) in enumerate(zip(self.kpt_coords, self.kpt_weights), 1):
            kpt_block.append([i, pretty_coords(c), f"{w:4.3f}"])

        lines.append(tabulate(kpt_block, tablefmt="plain"))
        lines.append("")
        lines.append("Band info near band edges")

        band_block = [["Index", "Kpoint index", "Energy", "Occupation",
                       "P-ratio", "Orbital"]]
        for oi in self.orbital_infos:
            orb_infos = np.array(oi).T

            # determine the band_idx where the occupation changes
            middle_idx = int(len(orb_infos) / 2)
            for band_idx in range(0, len(orb_infos) - 1):
                occu_diff = (orb_infos[band_idx][0].occupation
                             - orb_infos[band_idx + 1][0].occupation)
                if occu_diff > 0.1:
                    middle_idx = band_idx
                    break
            min_idx = max(middle_idx - 3, 0)
            max_idx = min(middle_idx + 3, len(orb_infos))

            for band_idx in range(min_idx, max_idx):
                for kpt_idx, orb_info in enumerate(orb_infos[band_idx], 1):
                    actual_band_idx = band_idx + self.lowest_band_index + 1
                    energy = f"{orb_info.energy:5.2f}"
                    occupation = f"{orb_info.occupation:4.1f}"
                    p_ratio = f"{orb_info.participation_ratio:4.1f}"
                    orbs = pretty_orbital(orb_info.orbitals)
                    band_block.append([actual_band_idx, kpt_idx, energy,
                                       occupation, p_ratio, orbs])
                band_block.append(["--"])
            band_block.append("")
        lines.append(tabulate(band_block, tablefmt="plain"))

        return "\n".join(lines)


def is_any_part_occ(x: OrbitalInfo, y: OrbitalInfo):
    return (min(x.occupation, (1.0 - x.occupation)) > 0.01
            or min(y.occupation, (1.0 - y.occupation)) > 0.01)


@dataclass
class LocalizedOrbital(MSONable):
    band_idx: int
    ave_energy: float
    occupation: float
    orbitals: Dict[str, List[float]]
    participation_ratio: Optional[float] = None
    radius: Optional[float] = None
    center: Optional[Coords] = None


@dataclass
class EdgeInfo(MSONable):
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
    vbm_info: EdgeInfo  # None for metals
    cbm_info: EdgeInfo

    def __str__(self):
        inner_table = [["", "Index", "Energy", "Occupation", "Orbitals", "K-point coords"]]
        inner_table.append(["VBM", self.vbm_info.band_idx, self.vbm_info.energy,
                            f"{self.vbm_info.occupation:5.2f}",
                            pretty_orbital(self.vbm_info.orbital_info.orbitals),
                            pretty_coords(self.vbm_info.kpt_coord),
                            ])
        inner_table.append(["CBM", self.cbm_info.band_idx, self.cbm_info.energy,
                            f"{self.cbm_info.occupation:5.2f}",
                            pretty_orbital(self.cbm_info.orbital_info.orbitals),
                            pretty_coords(self.cbm_info.kpt_coord),
                            ])
        return tabulate(inner_table, tablefmt="plain")


@dataclass
class BandEdgeState(MSONable):
    vbm_info: EdgeInfo
    cbm_info: EdgeInfo
    vbm_orbital_diff: float
    cbm_orbital_diff: float
    localized_orbitals: List[LocalizedOrbital]

    @property
    def is_shallow(self):
        return self.vbm_info.occupation < 0.7 or self.cbm_info.occupation > 0.3

    def __str__(self):
        lines = []
        inner_table = [["", "Index", "Energy", "P-ratio", "Occupation",
                        "Orbitals", "K-point coords"]]
        inner_table.append(["VBM", self.vbm_info.band_idx,
                            f"{self.vbm_info.energy:7.3f}",
                            f"{self.vbm_info.p_ratio:5.2f}",
                            f"{self.vbm_info.occupation:5.2f}",
                            pretty_orbital(self.vbm_info.orbital_info.orbitals),
                            pretty_coords(self.vbm_info.kpt_coord),
                            ])
        inner_table.append(["CBM", self.cbm_info.band_idx,
                            f"{self.cbm_info.energy:7.3f}",
                            f"{self.cbm_info.p_ratio:5.2f}",
                            f"{self.cbm_info.occupation:5.2f}",
                            pretty_orbital(self.cbm_info.orbital_info.orbitals),
                            pretty_coords(self.cbm_info.kpt_coord),
                            ])
        lines.append(tabulate(inner_table, tablefmt="plain"))

        lines.append("---")
        lines.append("Localized Orbital(s)")
        inner_table = [["Index", "Energy", "P-ratio", "Occupation", "Orbitals"]]
        for lo in self.localized_orbitals:
            inner_table.append([lo.band_idx,
                                f"{lo.ave_energy:7.3f}",
                                f"{lo.participation_ratio:5.2f}",
                                f"{lo.occupation:5.2f}",
                                pretty_orbital(lo.orbitals)])
        lines.append(tabulate(inner_table, tablefmt="plain"))
        return "\n".join(lines)


@dataclass
class BandEdgeStates(MSONable, ToJsonFileMixIn):
    states: List[BandEdgeState]  # by spin.

    @property
    def is_shallow(self):
        if any([i.is_shallow for i in self.states]):
            return IsShallow(True)
        return IsShallow(False)

    @property
    def band_indices(self):
        result = set()
        for state in self.states:
            result.add(state.vbm_info.band_idx)
            for lo in state.localized_orbitals:
                result.add(lo.band_idx)
            result.add(state.cbm_info.band_idx)
        return sorted(result)

    def __str__(self):
        lines = [" -- band-edge states info"]
        for spin, state in zip(["up", "down"], self.states):
            lines.append(f"Spin-{spin}")
            lines.append(state.__str__())
            lines.append("")

        return "\n".join(lines)


@dataclass
class IsShallow:
    is_shallow: bool

    def to_yaml(self, filename="is_shallow.yaml"):
        Path(filename).write_text(str(self.is_shallow))

    @classmethod
    def from_yaml(cls, filename="is_shallow.yaml"):
        first_letter = Path(filename).read_text()[0].upper()
        is_shallow = True if first_letter == "T" else False
        return cls(is_shallow=is_shallow)

