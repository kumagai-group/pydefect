# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from itertools import zip_longest
from typing import Dict, List

from pydefect.analyzer.defect_energy import CrossPoints
from tabulate import tabulate


@dataclass
class TransitionLevel:
    name: str
    charges: List[List[int]]  # [[2, 1], [1, 0]]
    energies: List[float]
    fermi_levels: List[float]


@dataclass
class TransitionLevels:
    transition_levels: List[TransitionLevel]
    cbm: float
    supercell_vbm: float
    supercell_cbm: float


def make_transition_levels(cross_point_dicts: Dict[str, CrossPoints],
                           cbm: float,
                           supercell_vbm: float,
                           supercell_cbm: float) -> TransitionLevels:
    transition_levels = []
    for name, cross_point in cross_point_dicts.items():
        charges = [list(i) for i in cross_point.charge_list[1:-1]]
        energies = cross_point.t_inner_cross_points[1]
        fermi_levels = cross_point.t_inner_cross_points[0]
        transition_levels.append(
            TransitionLevel(name, charges, energies, fermi_levels))

    return TransitionLevels(
        transition_levels, cbm, supercell_vbm, supercell_cbm)


def show_transition_levels(transition_levels: TransitionLevels):
    result = []
    for tl in transition_levels.transition_levels:
        for name, charge, energy, fermi in zip_longest(
                [tl.name], tl.charges, tl.energies, tl.fermi_levels,
                fillvalue=""):
            pretty_charges = f"{charge[0]} | {charge[1]}"
            result.append([name, pretty_charges, fermi, energy])
    headers = ["name", "charges", "Fermi level", "Formation energy"]
    return tabulate(result, headers=headers)