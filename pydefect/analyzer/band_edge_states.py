# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

from monty.json import MSONable
from vise.util.enum import ExtendedEnum

from pydefect.util.mix_in import ToJsonFileMixIn


@dataclass
class BandEdgeEigenvalues(MSONable, ToJsonFileMixIn):
    energies_and_occupations: List[List[List[float]]]  # [spin, k-idx, band-idx]
    kpt_coords: List[Tuple[float, float, float]]


@dataclass
class BandEdgeStates(MSONable, ToJsonFileMixIn):
    states: List["EdgeState"]  # by spin.
    method: str

    @property
    def is_shallow(self):
        return any([s.is_shallow for s in self.states])


class EdgeState(MSONable, ExtendedEnum):
    donor_phs = "Donor PHS"
    acceptor_phs = "Acceptor PHS"
    in_gap_state = "In-gap state"
    no_in_gap = "No in-gap state"
    unknown = "Unknown"

    @property
    def is_shallow(self):
        return self in [self.acceptor_phs, self.donor_phs]


@dataclass
class EdgeCharacters(MSONable, ToJsonFileMixIn):
    edge_characters: List["EdgeCharacter"]  # [by spin]


@dataclass
class EdgeCharacter(MSONable):
    """Code and its version dependent quantities. """
    hob_bottom_e: float  # hob bottom energy
    lub_top_e: float  # lub top energy
    vbm: Optional[float]
    cbm: Optional[float]
    # {"Mn": [0.01, ..], "O": [0.03, 0.5]},
    # where lists contain s, p, d, (f) orbital components.
    vbm_orbitals: Dict[str, List[float]]
    cbm_orbitals: Dict[str, List[float]]
    hob_participation_ratio: Optional[float] = None  # participation ratio
    lub_participation_ratio: Optional[float] = None
