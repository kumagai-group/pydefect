# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Dict, Tuple

import numpy as np
from monty.json import MSONable
from vise.util.enum import ExtendedEnum

from pydefect.util.mix_in import ToJsonFileMixIn


@dataclass
class BandEdgeEigenvalues(MSONable, ToJsonFileMixIn):
    energies: np.ndarray  # [spin, k-idx, band-idx]
    occupations: np.ndarray  # [spin, k-idx, band-idx]
    lowest_band_energy: float  # top of lowest band energy in this data
    highest_band_energy: float  # bottom of highest band energy in this data
    kpt_coords: List[Tuple[float, float, float]]
    kpt_weights: List[float]


@dataclass
class BandEdgeStates(MSONable, ToJsonFileMixIn):
    states: List["EdgeState"]  # by spin.
    method: str


class EdgeState(MSONable, ExtendedEnum):
    donor_phs = "Donor PHS"
    acceptor_phs = "Acceptor PHS"
    localized_state = "Localized state"
    no_in_gap = "No in-gap state"

    @property
    def is_shallow(self):
        return self in [self.acceptor_phs, self.donor_phs]


@dataclass
class EdgeCharacters(MSONable, ToJsonFileMixIn):
    # [by spin][hob, lub]
    edge_characters: List[Tuple["EdgeCharacter", "EdgeCharacter"]]


@dataclass
class EdgeCharacter(MSONable):
    """Code and its version dependent quantities. """
    participation_ratio: float
    energy: float
    # {"Mn": [0.01, ..], "O": [0.03, 0.5]},
    # where lists contain s, p, d, (f) orbital components.
    orbitals: Dict[str, List[float]]
