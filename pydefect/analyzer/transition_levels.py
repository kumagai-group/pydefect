# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from itertools import zip_longest
from typing import Dict, List

from monty.json import MSONable
from pydefect.analyzer.defect_energy import CrossPoints
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class TransitionLevel(MSONable):
    name: str
    charges: List[List[int]]  # [[2, 1], [1, 0]]
    energies: List[float]
    fermi_levels: List[float]


@dataclass
class TransitionLevels(MSONable, ToJsonFileMixIn):
    transition_levels: List[TransitionLevel]
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    def __str__(self):
        header = f"vbm: 0.00, cbm: {self.cbm:.2f}, " \
                 f"supercell vbm: {self.supercell_vbm:.2f}, " \
                 f"supercell_cbm: {self.supercell_cbm:.2f}\n"
        result = []
        for tl in self.transition_levels:
            if not tl.energies:
                continue
            for name, charge, energy, fermi in zip_longest(
                    [tl.name], tl.charges, tl.energies, tl.fermi_levels,
                    fillvalue=""):
                pretty_charges = f"{charge[0]} | {charge[1]}"
                result.append([name, pretty_charges, fermi, energy])
        headers = ["name", "charges", "Fermi level", "Formation energy"]
        floatfmt = ("", "", ".3f", ".3f", "")
        return header + tabulate(result, headers=headers, floatfmt=floatfmt)


def make_transition_levels(cross_point_dicts: Dict[str, CrossPoints],
                           cbm: float,
                           supercell_vbm: float,
                           supercell_cbm: float) -> TransitionLevels:
    transition_levels = []
    for name, cross_point in cross_point_dicts.items():
        charges = [list(i) for i in cross_point.charge_list[1:-1]]
        if cross_point.inner_cross_points:
            energies = cross_point.t_inner_cross_points[1]
            fermi_levels = cross_point.t_inner_cross_points[0]
        else:
            energies, fermi_levels = [], []
        transition_levels.append(
            TransitionLevel(name, charges, energies, fermi_levels))

    return TransitionLevels(
        transition_levels, cbm, supercell_vbm, supercell_cbm)

