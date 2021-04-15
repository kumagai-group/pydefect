# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy
from dataclasses import dataclass, InitVar
from itertools import groupby, combinations
from typing import List, Dict, Optional, Tuple

import dash_html_components as html
import numpy as np
from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen.core import Element
from scipy.spatial import HalfspaceIntersection
from vise.util.mix_in import ToYamlFileMixIn


@dataclass
class DefectEnergy(MSONable):
    formation_energy: float
    correction_energy: Dict[str, float]
    is_shallow: Optional[bool] = None

    @property
    def total_correction(self) -> float:
        return sum([v for v in self.correction_energy.values()])


@dataclass
class DefectEnergyInfo(MSONable, ToYamlFileMixIn):
    name: str
    charge: int
    atom_io: Dict[Element, int]
    energy: DefectEnergy

    def to_yaml(self) -> str:
        lines = [f"name: {self.name}",
                 f"charge: {self.charge}",
                 f"formation_energy: {self.energy.formation_energy}",
                 f"atom_io:"]
        for k, v in self.atom_io.items():
            lines.append(f"  {k}: {v}")
        lines.append(f"correction_energy:")
        for k, v in self.energy.correction_energy.items():
            lines.append(f"  {k}: {v}")
        is_shallow = "" if self.energy.is_shallow is None \
            else self.energy.is_shallow
        lines.append(f"is_shallow: {is_shallow}")
        return "\n".join(lines)

    @classmethod
    def from_yaml(cls, filename: str = "energy.yaml") -> "DefectEnergyInfo":
        d = loadfn(filename)
        d["atom_io"] = {Element(k): v for k, v in d["atom_io"].items()} \
            if d["atom_io"] else {}
        return cls(d.pop("name"), d.pop("charge"), d.pop("atom_io"),
                   DefectEnergy(**d))


@dataclass
class DefectEnergySummary:  # yaml
    title: str
    defect_energies: Dict[str, "DefectEnergies"]
    rel_chem_pot: Dict[str, Dict[Element, float]]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    def make_charge_and_energies(self, chem_pot, is_shallow, with_correction,
                                 band_edge):
        pass


@dataclass
class DefectEnergies:
    atom_io: Dict[Element, int]
    charges: List[int]
    energies: List[DefectEnergy]

    @property
    def corrected_energies(self):
        return [e.correction_energy for e in self.energies]

    def __str__(self):
        lines = []
        for c, e, cor in zip(self.charges, self.energies, self.corrections):
            lines.append(f"{self.name:>10} {c:>4} {e:12.4f} {cor:12.4f}")
        return "\n".join(lines)


@dataclass
class ChargesAndEnergies:
    charges: Dict[str, List[int]]
    energies: Dict[str, List[float]]
    ef_min: InitVar[float]
    ef_max: InitVar[float]
    base_e: float = 0.0

    def __post_init__(self, ef_min: float, ef_max: float):
        large_minus_number = -1e4
        half_spaces = []
        for charge, corr_energy in zip(self.charges, self.energies):
            half_spaces.append([-charge, 1, -corr_energy])

        half_spaces.append([-1, 0, ef_min])
        half_spaces.append([1, 0, -ef_max])
        half_spaces.append([0, -1, large_minus_number])

        feasible_point = np.array([(ef_min + ef_max) / 2, -1e3])

        hs = HalfspaceIntersection(np.array(half_spaces), feasible_point)
        boundary_points = []
        inner_cross_points = []
        for intersection in hs.intersections:
            x, y = np.round(intersection, 8)
            if ef_min + 0.001 < x < ef_max - 0.001:
                inner_cross_points.append([x - self.base_e, y])
            elif y > large_minus_number + 1:
                boundary_points.append([x - self.base_e, y])

        self.cross_points = CrossPoints(inner_cross_points, boundary_points)

    @property
    def stable_charges(self):
        return set(self.cross_points.charges)

    @property
    def transition_levels(self) -> Dict[Tuple[int, int], float]:
        result = {}
        for (c1, e1), (c2, e2) in combinations(
                zip(self.charges, self.energies), 2):
            result[(c1, c2)] = - (e1 - e2) / (c1 - c2) - self.base_e
        return result

    @property
    def pinning_level(self) -> Tuple[Tuple[float, Optional[int]],
                                     Tuple[float, Optional[int]]]:
        """
        :param base_e: Reference to show the pinning level such as VBM.
        :return: ((Lower pinning, its charge), (Upper pinning, its charge))
        """
        lower_pinning, upper_pinning = float("-inf"), float("inf")
        lower_charge, upper_charge = None, None
        for charge, corr_energy in zip(self.charges, self.energies):
            if charge == 0:
                continue
            pinning = - corr_energy / charge
            if charge > 0 and pinning > lower_pinning:
                lower_pinning, lower_charge = pinning, charge
            elif pinning < upper_pinning:
                upper_pinning, upper_charge = pinning, charge
        return ((lower_pinning - self.base_e, lower_charge),
                (upper_pinning - self.base_e, upper_charge))

    def energy_at_ef(self, ef: float) -> Tuple[float, int]:
        """
        :return: (Lowest energy, its charge)
        """
        result_e, result_charge = float("inf"), None
        for charge, corr_energy in zip(self.charges, self.energies):
            energy = corr_energy + charge * ef
            if energy < result_e:
                result_e, result_charge = energy, charge
        return result_e, result_charge


@dataclass
class CrossPoints:
    inner_cross_points: List[List[float]]
    boundary_points: List[List[float]]

    @property
    def all_sorted_points(self):
        return sorted(self.boundary_points + self.inner_cross_points,
                      key=lambda v: v[0])

    @property
    def t_all_sorted_points(self):
        return np.transpose(np.array(self.all_sorted_points)).tolist()

    @property
    def t_inner_cross_points(self):
        return np.transpose(np.array(self.inner_cross_points)).tolist()

    @property
    def t_boundary_points(self):
        return np.transpose(np.array(self.boundary_points)).tolist()

    @property
    def charges(self) -> List[int]:
        result = []
        for i, j in zip(self.all_sorted_points[:-1], self.all_sorted_points[1:]):
            dx = j[0] - i[0]
            dy = j[1] - i[1]
            result.append(int(round(dy / dx)))
        return result

    @property
    def charge_list(self):
        charges = [None] + self.charges + [None]
        return list(zip(charges[:-1], charges[1:]))

    @property
    def annotated_charge_positions(self):
        result = {}
        for ((x1, y1), (x2, y2)), charge \
                in zip(zip(self.all_sorted_points[:-1],
                           self.all_sorted_points[1:]),
                       self.charges):
            result[charge] = [(x1 + x2) / 2, (y1 + y2) / 2]
        return result

    def __str__(self):
        lines = []
        for point in self.all_sorted_points:
            lines.append(f"{point[0]:12.4f} {point[1]:12.4f}")
        return "\n".join(lines)


def make_defect_energies(single_energies: List[DefectEnergyInfo],
                         abs_chem_pot:  Dict[Element, float],
                         allow_shallow: bool,
                         ) -> List[DefectEnergy]:
    if allow_shallow is False:
        single_energies = [e for e in single_energies if e.is_shallow is False]
    sorted_energies = sorted(single_energies, key=lambda x: x.name)
    result = []
    for _, grouped_es in groupby(sorted_energies, lambda x: x.name):
        grouped_es = list(grouped_es)
        name = grouped_es[0].name
        charges = [e.charge for e in grouped_es]
        energies = [e.formation_energy_wo_corr(abs_chem_pot) for e in grouped_es]
        corrections = [e.total_correction for e in grouped_es]
        result.append(DefectEnergy(name, charges, energies, corrections))
    return result


def remove_digits(name):
    return ''.join([i for i in name if not i.isdigit()])


def only_digits(name):
    return ''.join([i for i in name if i.isdigit()])


elements = [str(e) for e in Element]


def defect_mpl_name(name):
    in_name, out_name = name.split("_")
    if in_name in elements:
        in_name = "{\\rm " + in_name + "}"
    elif in_name == "Va":
        in_name = "V"

    r_out_name = remove_digits(out_name)
    if r_out_name in elements:
        out_name = "{{\\rm " + r_out_name + "}" + only_digits(out_name) + "}"
    else:
        out_name = "{" + out_name + "}"

    return f"${in_name}_{out_name}$"


def defect_plotly_name(name):
    in_name, out_name = name.split("_")
    if in_name == "Va":
        in_name = "<i>V</i>"

    out_name = f"<sub>{out_name}</sub>"
    return f"{in_name}{out_name}"


def defect_plotly_full_name(fullname):
    in_name, out_name, charge = fullname.split("_")
    if in_name == "Va":
        in_name = "<i>V</i>"
    return f"{in_name}<sub>{out_name}</sub><sup>{charge}</sup>"


def defect_html_title_name(fullname):
    x = fullname.split("_")
    if len(x) == 2:
        in_name, out_name = x
    elif len(x) == 3:
        in_name, out_name, charge = x
    else:
        raise ValueError

    if in_name == "Va":
        in_name = html.I("V")
    else:
        in_name = html.Span(in_name)

    result = [in_name, html.Sub(out_name)]
    if len(x) == 3:
        result.append(html.Sup(charge))
    return result


def sanitize_defect_energies_for_plot(defect_energies: List[DefectEnergy],
                                      for_plotly: bool = False):
    result = []
    out_names = [e.name.split("_")[1] for e in defect_energies]

    for e in defect_energies:
        ee = deepcopy(e)
        in_name, out_name = e.name.split("_")
        r_out_name = remove_digits(out_name)
        out_name = r_out_name if f"{r_out_name}2" not in out_names else out_name
        if for_plotly:
            ee.name = defect_plotly_name("_".join([in_name, out_name]))
        else:
            ee.name = defect_mpl_name("_".join([in_name, out_name]))
        result.append(ee)

    return result


def slide_energy(defect_energies: List[DefectEnergy], base_level: float):
    result = []
    for e in defect_energies:
        ee = deepcopy(e)
        ee.energies = [e + base_level * c for e, c in zip(e.energies, e.charges)]
        result.append(ee)
    return result


def reservoir_energy(atom_io: Dict[Element, int],
                     abs_chem_pot: Dict[Element, float]) -> float:
    return sum([-diff * abs_chem_pot[elem] for elem, diff in atom_io.items()])