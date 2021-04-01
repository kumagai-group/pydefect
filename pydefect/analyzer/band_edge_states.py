# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections.abc import Sequence
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np

from monty.json import MSONable
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class BandEdgeEigenvalues(MSONable, ToJsonFileMixIn):
    # [spin, k-idx, band-idx] = energy, occupation
    energies_and_occupations: List[List[List[List[float]]]]
    kpt_coords: List[Tuple[float, float, float]]
    lowest_band_index: int


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
    kpt_coords: List[Tuple[float, float, float]]
    kpt_weights: List[float]
    lowest_band_index: int

    @property
    def energies_and_occupations(self):
        result = deepcopy(self.orbital_infos)
        for i, x in enumerate(self.orbital_infos):
            for j, y in enumerate(x):
                for k, z in enumerate(y):
                    result[i][j][k] = [z.energy, z.occupation]
        return result

    def __repr__(self):
        lines = ["k-points"]
        x = [["idx", "coords", "weight"]]
        for i, (c, w) in enumerate(zip(self.kpt_coords, self.kpt_weights), 1):
            x.append([i, c, w])
        lines.append(tabulate(x))
        x = [["band_idx", "k_idx", "energy", "occupation", "p_ratio", "orbital"]]
        for oi in self.orbital_infos:
            a = np.array(oi).T
            for i, b in enumerate(a, 1):
                for j, c in enumerate(b, 1):
                    x.append([i+self.lowest_band_index, j, c.energy, round(c.occupation, 3), round(c.participation_ratio, 2), c.orbitals])
            x.append("")
        lines.append(tabulate(x))

        return "\n".join(lines)


@dataclass
class LocalizedOrbital(MSONable):
    band_idx: int
    ave_energy: float
    occupation: float
    orbitals: Dict[str, List[float]]
#    radius: Optional[float] = None
#    center: Optional[Tuple[float, float, float]] = None


@dataclass
class EdgeInfo(MSONable):
    band_idx: int
    kpt_coord: Tuple[float, float, float]
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


@dataclass
class PerfectBandEdgeState(MSONable, ToJsonFileMixIn):
    vbm_info: EdgeInfo  # None for metals
    cbm_info: EdgeInfo


@dataclass
class BandEdgeState(MSONable):
    vbm_info: EdgeInfo
    cbm_info: EdgeInfo
    vbm_orbital_diff: float
    cbm_orbital_diff: float
    localized_orbitals: List[LocalizedOrbital]

    @property
    def is_shallow(self):
        return self.vbm_info.occupation < 0.5 or self.cbm_info.occupation > 0.5


@dataclass
class BandEdgeStates(MSONable, ToJsonFileMixIn):
    states: List[BandEdgeState]  # by spin.

    @property
    def is_shallow(self):
        if any([i.is_shallow for i in self.states]):
            return IsShallow(True)
        return IsShallow(False)


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

