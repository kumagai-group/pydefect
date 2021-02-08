# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import math
from dataclasses import dataclass
from typing import Optional, List, Tuple

import numpy as np
from monty.json import MSONable
from numpy.linalg import inv
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator, SiteDiff
from pymatgen import Structure, Lattice
from pymatgen.symmetry.groups import SpaceGroup
from tabulate import tabulate
from vise.util.enum import ExtendedEnum
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


def fold_coords(structure: Structure, center: List[float]):
    for site in structure:
        _, image = site.distance_and_image_from_frac_coords(center)
        site.frac_coords -= image


def calc_drift(perfect: Structure, defect: Structure, center: List[float],
               d_to_p: List[int]):
    distances = []
    for site in defect:
        distances.append(site.distance_and_image_from_frac_coords(center)[0])
    anchor_atom_idx = int(np.argmax(distances))
    p_anchor_atom_idx = d_to_p[anchor_atom_idx]
    distance = defect[anchor_atom_idx].distance(perfect[p_anchor_atom_idx])
    drift_vector = (defect[anchor_atom_idx].frac_coords
                    - perfect[p_anchor_atom_idx].frac_coords)
    return anchor_atom_idx, distance, drift_vector


def calc_displacements(perfect: Structure,
                       defect: Structure,
                       center: List[float],
                       d_to_p: List[int],
                       target_idxs: List[int] = None):
    result = []
    lattice: Lattice = defect.lattice
    for d, p in enumerate(d_to_p):
        if target_idxs and d not in target_idxs:
            continue
        if p is None:
            result.append(None)
        else:
            p_elem = str(perfect[p].specie)
            d_elem = str(defect[d].specie)
            if p_elem != d_elem:
                result.append(None)
                continue
            orig_pos = tuple(perfect[p].frac_coords)
            final_pos = tuple(defect[d].frac_coords)

            initial_pos_vec = \
                lattice.get_cartesian_coords(perfect.frac_coords[p] - center)
            disp_dist, t = lattice.get_distance_and_image(
                defect.frac_coords[d], perfect.frac_coords[p])
            displacement = lattice.get_cartesian_coords(
                defect.frac_coords[d] - perfect.frac_coords[p] + t)

            inner_prod = sum(initial_pos_vec * displacement)
            ini_dist = np.linalg.norm(initial_pos_vec)
            cos = round(inner_prod / (ini_dist * disp_dist), 10)
            angle = float(round(180 * (1 - np.arccos(cos) / np.pi), 1))
            if math.isnan(angle):
                angle = None

            result.append(Displacement(specie=p_elem,
                                       original_pos=orig_pos,
                                       final_pos=final_pos,
                                       distance_from_defect=ini_dist,
                                       displace_distance=disp_dist,
                                       angle=angle))
    return result


def remove_dot(x):
    return "".join([s for s in x if s != "."])


def elem_indices_coords(structure: Structure, indices):
    result = []
    for i in indices:
        site = structure[i]
        result.append((str(site.specie), i, tuple(site.frac_coords)))
    return result

# fold_coords(perfect, center)
# fold_coords(initial, center)
# fold_coords(final, center)


def make_defect_structure_info(perfect: Structure,
                               initial: Structure,
                               final: Structure,
                               dist_tol: float,
                               symprec: float):
    final = final.copy()
    initial_symmetrizer = StructureSymmetrizer(initial, symprec)
    final_symmetrizer = StructureSymmetrizer(final, symprec)

    assert perfect.lattice == initial.lattice == final.lattice
    comp_w_perf = DefectStructureComparator(final, perfect, dist_tol)
    comp_w_init = DefectStructureComparator(final, initial, dist_tol)
    d_to_p = comp_w_init.d_to_p

    center = comp_w_perf.defect_center_coord
    anchor_idx, drift_d, drift_vec = calc_drift(perfect, final, center, d_to_p)
    center -= drift_vec
    for site in final:
        site.frac_coords -= drift_vec

    displacements = calc_displacements(perfect, final, center, d_to_p)

    return DefectStructureInfo(
        initial_site_sym=unique_point_group(initial_symmetrizer.point_group),
        final_site_sym=unique_point_group(final_symmetrizer.point_group),
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
class DefectStructureInfo(MSONable):
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

    def __repr__(self):
        lines = [f"Site symmetry: {self.initial_site_sym} "
                 f"-> {self.final_site_sym} ({self.symm_relation})",
                 f"Is same configuration: {self.same_config_from_init}",
                 f"Drift distance: {self.drift_dist:5.3f}"]
        center_coords = [round(c, 3) for c in self.center]
        lines.append(f"Defect center: {center_coords}")

        lines.append("Removed atoms:")
        x = []
        for i, (elem, coords) in self.site_diff.removed.items():
            x.append([i, elem] + [round(fc, 2) for fc in coords])
        lines.append(tabulate(x))

        lines.append("Added atoms:")
        x = []
        for i, (elem, coords) in self.site_diff.inserted.items():
            x.append([i, elem] + [round(fc, 2) for fc in coords])
        lines.append(tabulate(x))

        if self.same_config_from_init is False:
            lines.append("Initially removed atoms:")
            x = []
            for i, (elem, coords) in self.site_diff_from_initial.removed.items():
                x.append([i, elem] + [round(fc, 2) for fc in coords])
            lines.append(tabulate(x))

            lines.append("Initially added atoms:")
            x = []
            for i, (elem, coords) in self.site_diff_from_initial.inserted.items():
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


    # def to(self, dir_name: Path = Path(".")):
    #     fd_local = Structure.from_sites(
    #         [self.final_structure[s] for s in self.neighboring_atom_indices])
    #     for v in self.final_vacancies:
    #         fd_local.append("X", coords=self.perfect_structure[v].frac_coords)
    #     fd_local.to(filename=dir_name / "POSCAR_final_local")
    #
    #     neighbor_indices_in_p = [self.fd_to_p[i]
    #                              for i in self.neighboring_atom_indices
    #                              if self.fd_to_p[i] is not None]
    #     neighbor_indices_in_p.extend(self.final_vacancies)
    #     p_local = Structure.from_sites(
    #         [self.perfect_structure[s] for s in neighbor_indices_in_p])
    #     p_local.to(filename=dir_name / "POSCAR_perfect_local")
    #
