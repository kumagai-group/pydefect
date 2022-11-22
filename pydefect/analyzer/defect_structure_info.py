# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import Optional, List, Tuple

from monty.json import MSONable
from pydefect.analyzer.defect_structure_comparator import \
    SiteDiff, SiteInfo
from pydefect.defaults import defaults
from pydefect.util.coords import pretty_coords
from pymatgen.core import Structure
from pymatgen.symmetry.groups import SpaceGroup
from tabulate import tabulate
from vise.util.enum import ExtendedEnum
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn

logger = get_logger(__name__)


class DefectType(MSONable, ExtendedEnum):
    vacancy = "vacancy"
    interstitial = "interstitial"
    substituted = "substituted"
    vacancy_split = "vacancy_split"
    interstitial_split = "interstitial_split"
    unknown = "unknown"


def judge_defect_type(site_diff: SiteDiff):
    if site_diff.is_vacancy:
        return DefectType.vacancy
    elif site_diff.is_interstitial:
        return DefectType.interstitial
    elif site_diff.is_substituted:
        return DefectType.substituted

    elements_involved = set()
    for _, elem, _ in site_diff.removed + site_diff.inserted:
        elements_involved.add(elem)

    if len(elements_involved) == 1 and not site_diff.removed_by_sub:
        if len(site_diff.removed) - len(site_diff.inserted) == 1:
            return DefectType.vacancy_split
        elif len(site_diff.removed) - len(site_diff.inserted) == -1:
            return DefectType.interstitial_split

    return DefectType.unknown


def remove_dot(x):
    return "".join([s for s in x if s != "."])


def unique_point_group(pg):
    result = remove_dot(pg)
    if result == "2mm" or result == "m2m":
        return "mm2"
    if result == "-4m2":
        return "-42m"
    if result == "m3":
        return "m-3"
    return result


class SymmRelation(MSONable, ExtendedEnum):
    same = "same"
    subgroup = "subgroup"
    supergroup = "supergroup"
    another = "another"


@dataclass
class Displacement(MSONable):
    specie: str
    original_pos: Tuple[float, float, float]
    final_pos: Tuple[float, float, float]
    distance_from_defect: float
    disp_vector: Tuple[float, float, float]
    displace_distance: float
    angle: Optional[float]


def symmetry_relation(initial_point_group, final_point_group):
    """ Check the point group symmetry relation using the space group relation
    implemented in pymatgen.
    """
    if initial_point_group in ["3m", "-3m"]:
        initial_point_group += "1"
    if final_point_group in ["3m", "-3m"]:
        final_point_group += "1"

    initial = SpaceGroup(f"P{initial_point_group}")
    final = SpaceGroup(f"P{final_point_group}")
    if initial == final:
        return SymmRelation.same
    elif final.is_subgroup(initial):
        return SymmRelation.subgroup
    elif final.is_supergroup(initial):
        return SymmRelation.supergroup
    else:
        return SymmRelation.another


@dataclass
class DefectStructureInfo(MSONable, ToJsonFileMixIn):
    shifted_final_structure: Structure
    initial_site_sym: str
    final_site_sym: str
    site_diff: SiteDiff
    site_diff_from_initial: SiteDiff
    symprec: float
    dist_tol: float
    anchor_atom_idx: int
    neighbor_atom_indices: List[int]
    neighbor_cutoff_factor: float
    drift_vector: Tuple[float, float, float]
    drift_dist: float
    center: Tuple[float, float, float]
    displacements: List[Displacement]

    @property
    def symm_relation(self):
        return symmetry_relation(self.initial_site_sym, self.final_site_sym)

    @property
    def same_config_from_init(self):
        return self.site_diff_from_initial.is_no_diff

    @property
    def defect_type(self):
        return judge_defect_type(self.site_diff)

    def __str__(self):
        sym_transition = f"{self.initial_site_sym} " \
                         f"-> {self.final_site_sym} ({self.symm_relation})"
        center_coords = pretty_coords(self.center)
        drift_dist = f"{self.drift_dist:5.3f}" if self.drift_dist else "N.A."
        lines = [" -- defect structure info",
                 f"Defect type: {self.defect_type}",
                 f"Site symmetry: {sym_transition}",
                 f"Has same configuration from initial structure: "
                 f"{self.same_config_from_init}",
                 f"Drift distance: {drift_dist}",
                 f"Defect center: {center_coords}"]

        def _site_info(header: str, site_info: List[SiteInfo]):
            if not site_info:
                return []
            _table = tabulate([[idx, elem, pretty_coords(coords)]
                              for idx, elem, coords in site_info],
                              tablefmt="plain")
            return [header, _table, ""]

        lines.extend(_site_info("Removed atoms:", self.site_diff.removed))
        lines.extend(_site_info("Added atoms:", self.site_diff.inserted))

        if self.same_config_from_init is False:
            lines.extend(_site_info("Removed atoms from initial structure:",
                                    self.site_diff_from_initial.removed))
            lines.extend(_site_info("Inserted atoms to initial structure:",
                                    self.site_diff_from_initial.inserted))

        min_d = min([d.distance_from_defect for d in self.displacements if d])
        cutoff = min_d * self.neighbor_cutoff_factor
        lines.append(f"Neighbor max distance {cutoff:5.3f}")

        lines.append("Displacements")
        idxs = [[i, d] for i, d in enumerate(self.displacements)
                if d is not None]
        table = [["Elem", "Dist", "Displace", "Angle", "Index",
                  "Initial site", "", "Final site", "Neighbor"]]
        print(self.neighbor_atom_indices)
        for final_idx, d in sorted(idxs,
                                   key=lambda y: y[1].distance_from_defect):
            is_neighbor = "T" if final_idx in self.neighbor_atom_indices else ""
            if d.distance_from_defect > defaults.show_structure_cutoff:
                break
            i_pos = pretty_coords(d.original_pos)
            f_pos = pretty_coords(d.final_pos)
            angle = int(round(d.angle, -1)) if d.angle else ""
            table.append([d.specie, round(d.distance_from_defect, 2),
                          round(d.displace_distance, 2), angle, final_idx,
                          i_pos, "->", f_pos, is_neighbor])
        lines.append(tabulate(table, tablefmt="plain"))

        return "\n".join(lines)


