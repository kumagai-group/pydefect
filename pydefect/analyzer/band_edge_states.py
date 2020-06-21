# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections.abc import Sequence
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

from monty.json import MSONable
from pydefect.util.mix_in import ToJsonFileMixIn
from vise.util.enum import ExtendedEnum


@dataclass
class BandEdgeEigenvalues(MSONable, ToJsonFileMixIn):
    energies_and_occupations: List[List[List[List[float]]]]  # [spin, k-idx, band-idx] = energy, occupation
    kpt_coords: List[Tuple[float, float, float]]
    lowest_band_index: int


class EdgeState(MSONable, ExtendedEnum):
    donor_phs = "Donor PHS"
    acceptor_phs = "Acceptor PHS"
    in_gap_state = "In-gap state"
    no_in_gap = "No in-gap state"
    unknown = "Unknown"

    @property
    def is_shallow(self):
        return self in [self.acceptor_phs, self.donor_phs]

    def as_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "value": self.value}


@dataclass
class BandEdgeStates(MSONable, ToJsonFileMixIn):
    states: List[EdgeState]  # by spin.

    @property
    def is_shallow(self):
        return any([s.is_shallow for s in self.states])


class EdgeCharacters(Sequence, MSONable, ToJsonFileMixIn):
    def __init__(self, edge_characters: List["EdgeCharacter"]):
        self.list = edge_characters  # [by spin]

    def __getitem__(self, item):
        return self.list[item]

    def __len__(self):
        return len(self.list)

    def as_dict(self):
        return {"edge_characters": [ec.as_dict() for ec in self.list]}

    @classmethod
    def from_dict(cls, d):
        return cls([EdgeCharacter.from_dict(dd) for dd in d["edge_characters"]])

    def __str__(self):
        strings = []
        for spin, edge_char in zip(["up", "down"], self.list):
            strings.append(f"""{spin:>4}: {edge_char}""")
        return "\n".join(strings)


@dataclass
class EdgeCharacter(MSONable):
    """Code and its version dependent quantities. """
    hob_bottom_e: float  # Highest-occupied band
    lub_top_e: float  # Lowest-unoccupied band
    vbm: Optional[float]  # None for metals
    cbm: Optional[float]
    # {"Mn": [0.01, ..], "O": [0.03, 0.5]},
    # where lists contain s, p, d, (f) orbital components.
    vbm_orbitals: Dict[str, List[float]]
    cbm_orbitals: Dict[str, List[float]]
    vbm_participation_ratio: Optional[float] = None  # participation ratio
    cbm_participation_ratio: Optional[float] = None
