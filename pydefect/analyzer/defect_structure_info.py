# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import math
import warnings
from dataclasses import dataclass
from typing import Optional, List, Tuple

import numpy as np
from monty.json import MSONable
from numpy.linalg import inv
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator, SiteDiff
from pymatgen.core import Structure, Lattice, PeriodicSite
from pymatgen.symmetry.groups import SpaceGroup
from tabulate import tabulate
from vise.util.enum import ExtendedEnum
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


def folded_coords(site: PeriodicSite, center: List[float]):
    _, image = site.distance_and_image_from_frac_coords(center)
    return tuple(site.frac_coords - image)


def fold_coords_in_structure(structure: Structure, center: List[float]):
    for site in structure:
        _, image = site.distance_and_image_from_frac_coords(center)
        site.frac_coords -= image


def calc_drift(perfect: Structure,
               defect: Structure,
               center: List[float],
               d_to_p: List[int]):
    distances = []
    for site in defect:
        distances.append(site.distance_and_image_from_frac_coords(center)[0])
    anchor_atom_idx = int(np.argmax(distances))
    p_anchor_atom_idx = d_to_p[anchor_atom_idx]
    distance, image = defect[anchor_atom_idx].distance_and_image_from_frac_coords(perfect[p_anchor_atom_idx].frac_coords)
    drift_vector = (defect[anchor_atom_idx].frac_coords
                    - perfect[p_anchor_atom_idx].frac_coords - image)
    return anchor_atom_idx, distance, drift_vector


def calc_displacements(perfect: Structure,
                       defect: Structure,
                       center: List[float],
                       d_to_p: List[int]):
    result = []
    lattice: Lattice = defect.lattice
    for d, p in enumerate(d_to_p):
        if p is None:
            result.append(None)
        else:
            p_elem = str(perfect[p].specie)
            d_elem = str(defect[d].specie)
            if p_elem != d_elem:
                result.append(None)
                continue
            initial_pos = folded_coords(perfect[p], center)
            final_pos = folded_coords(defect[d], initial_pos)

            initial_pos_vec = lattice.get_cartesian_coords(np.array(initial_pos) - center)
            disp_dist, t = lattice.get_distance_and_image(
                defect.frac_coords[d], perfect.frac_coords[p])
            disp_vec = lattice.get_cartesian_coords(
                defect.frac_coords[d] - perfect.frac_coords[p] - t)

            inner_prod = sum(initial_pos_vec * disp_vec)
            ini_dist = np.linalg.norm(initial_pos_vec)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                # ignore "RuntimeWarning: invalid value encountered in double_scalars"
                cos = round(inner_prod / (ini_dist * disp_dist), 10)
            angle = float(round(180 * (1 - np.arccos(cos) / np.pi), 1))
            if math.isnan(angle):
                angle = None

            result.append(Displacement(specie=p_elem,
                                       original_pos=initial_pos,
                                       final_pos=final_pos,
                                       distance_from_defect=ini_dist,
                                       disp_vector=tuple(disp_vec),
                                       displace_distance=disp_dist,
                                       angle=angle))
    return result


def remove_dot(x):
    return "".join([s for s in x if s != "."])


def elem_indices_coords(structure: Structure,
                        indices: List[int]
                        ) -> List[Tuple[str, int, Tuple[float, float, float]]]:
    result = []
    for i in indices:
        site = structure[i]
        result.append((str(site.specie), i, tuple(site.frac_coords)))
    return result


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


def make_defect_structure_info(perfect: Structure,
                               initial: Structure,
                               final: Structure,
                               dist_tol: float,
                               symprec: float = None,
                               init_site_sym: str = None,
                               final_site_sym: str = None
                               ) -> "DefectStructureInfo":
    final = final.copy()
    if symprec:
        initial_symmetrizer = StructureSymmetrizer(initial, symprec)
        final_symmetrizer = StructureSymmetrizer(final, symprec)
        init_site_sym = initial_symmetrizer.point_group
        final_site_sym = final_symmetrizer.point_group
    else:
        if not (init_site_sym and final_site_sym):
            raise AssertionError("Need to set both initial and final site "
                                 "symmetry, when symprec is not set.")

    assert perfect.lattice == initial.lattice == final.lattice
    comp_w_perf = DefectStructureComparator(final, perfect, dist_tol)
    comp_w_init = DefectStructureComparator(final, initial, dist_tol)
    d_to_p = comp_w_init.d_to_p

    center = comp_w_perf.defect_center_coord
    anchor_idx, drift_d, drift_vec = calc_drift(initial, final, center, d_to_p)
    center -= drift_vec
    for site in final:
        site.frac_coords -= drift_vec

    displacements = calc_displacements(initial, final, center, d_to_p)

    return DefectStructureInfo(
        initial_site_sym=unique_point_group(remove_dot(init_site_sym)),
        final_site_sym=unique_point_group(final_site_sym),
        site_diff=comp_w_perf.make_site_diff(),
        site_diff_from_initial=comp_w_init.make_site_diff(),
        symprec=symprec,
        dist_tol=dist_tol,
        anchor_atom_idx=anchor_idx,
        drift_vector=tuple(drift_vec),
        drift_dist=drift_d,
        center=tuple(center),
        displacements=displacements)


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


def unique_point_group(pg):
    if pg == "2mm" or pg == "m2m":
        return "mm2"
    if pg == "-4m2":
        return "-42m"
    if pg == "m3":
        return "m-3"
    return pg


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
    initial_site_sym: str
    final_site_sym: str
    site_diff: SiteDiff
    site_diff_from_initial: SiteDiff
    symprec: float
    dist_tol: float
    anchor_atom_idx: int
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

    def __repr__(self):
        lines = ["",
                 f"Defect type: {self.defect_type}",
                 f"Site symmetry: {self.initial_site_sym} "
                 f"-> {self.final_site_sym} ({self.symm_relation})",
                 f"Is same configuration: {self.same_config_from_init}",
                 f"Drift distance: {self.drift_dist:5.3f}"]
        center_coords = [round(c, 3) for c in self.center]
        lines.append(f"Defect center: {center_coords}")

        lines.append("Removed atoms:")
        x = []
        for i, elem, coords in self.site_diff.removed:
            x.append([i, elem] + [round(fc, 2) for fc in coords])
        lines.append(tabulate(x))

        lines.append("Added atoms:")
        x = []
        for i, elem, coords in self.site_diff.inserted:
            x.append([i, elem] + [round(fc, 2) for fc in coords])
        lines.append(tabulate(x))

        if self.same_config_from_init is False:
            lines.append("Initially removed atoms:")
            x = []
            for i, elem, coords in self.site_diff_from_initial.removed:
                x.append([i, elem] + [round(fc, 2) for fc in coords])
            lines.append(tabulate(x))

            lines.append("Initially added atoms:")
            x = []
            for i, elem, coords in self.site_diff_from_initial.inserted:
                x.append([i, elem] + [round(fc, 2) for fc in coords])
            lines.append(tabulate(x))

        lines.append("Displacements")
        idxs = [[i, d] for i, d in enumerate(self.displacements) if d is not None]
        x = []
        for final_idx, d in sorted(idxs, key=lambda y: y[1].distance_from_defect):
            i_pos = [round(c, 3) for c in d.original_pos]
            f_pos = [round(c, 3) for c in d.final_pos]
            angle = int(round(d.angle, -1)) if d.angle else ""
            x.append([d.specie, round(d.distance_from_defect, 2),
                      round(d.displace_distance, 2), angle, final_idx] + i_pos
                     + ["->"] + f_pos)
        lines.append(tabulate(x))

        return "\n".join(lines)


