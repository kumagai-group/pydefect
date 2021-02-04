# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
from monty.json import MSONable
from numpy.linalg import inv
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pydefect.input_maker.defect_entry import DefectEntry
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
    anchor_atom_idx = np.argmax(distances)
    p_anchor_atom_idx = d_to_p[anchor_atom_idx]

    distance = defect[anchor_atom_idx].distance(perfect[p_anchor_atom_idx])
    vector_in_frac = defect[anchor_atom_idx].frac_coords - perfect[p_anchor_atom_idx].frac_coords
    return distance, vector_in_frac


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
            elem = str(perfect[p].specie)
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

            result.append(Displacement(specie=elem,
                                       original_pos=orig_pos,
                                       final_pos=final_pos,
                                       distance_from_defect=ini_dist,
                                       displace_distance=disp_dist,
                                       angle=angle))
    return result


def remove_dot(x):
    return "".join([s for s in x if s != "."])


def make_defect_structure_info(perfect_calc_results: CalcResults,
                               defect_entry: DefectEntry,
                               calc_results: CalcResults,
                               neighbor_cutoff_factor: float = 1.5,
                               correct_drift: bool = False):
    perfect = perfect_calc_results.structure.copy()
    initial = defect_entry.structure.copy()
    final = calc_results.structure.copy()
    assert perfect.lattice == initial.lattice == final.lattice

    comp_i = DefectStructureComparator(initial, perfect)
    comp_f = DefectStructureComparator(final, perfect)
    neighbors = comp_f.neighboring_atom_indices(neighbor_cutoff_factor)
    center = comp_f.defect_center_coord

    drift_dist, drift_vec = calc_drift(perfect, final, center, comp_f.d_to_p)

    if correct_drift:
        center -= drift_vec
        for site in final:
            site.frac_coords -= drift_vec

    fold_coords(perfect, center)
    fold_coords(initial, center)
    fold_coords(final, center)

    displacements = calc_displacements(
        perfect, final, center, comp_f.d_to_p, neighbors)

    return DefectStructureInfo(initial_point_group=unique_poing_group(remove_dot(defect_entry.site_symmetry)),
                               final_point_group=unique_poing_group(calc_results.site_symmetry),
                               initial_structure=initial,
                               final_structure=final,
                               perfect_structure=perfect,
                               fd_to_p=comp_f.d_to_p,
                               drift_dist=drift_dist,
                               initial_vacancies=comp_i.vacant_indices,
                               initial_interstitials=comp_i.inserted_indices,
                               final_vacancies=comp_f.vacant_indices,
                               final_interstitials=comp_f.inserted_indices,
                               defect_center_coord=tuple(center),
                               displacements=displacements,
                               neighboring_atom_indices=neighbors)


class DefectType(MSONable, ExtendedEnum):
    vacancy = "vacancy"
    interstitial = "interstitial"
    substituted = "substituted"
    complex = "complex"


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


def unique_poing_group(pg):
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
    initial_point_group: str
    final_point_group: str
    initial_structure: Structure  # frac coords are shifted.
    final_structure: Structure  # frac coords are shifted.
    perfect_structure: Structure
    fd_to_p: List[int]
#    fd_to_id: List[int]
    drift_dist: float
#    initial_defect_type: DefectType
    initial_vacancies: List[int]
    initial_interstitials: List[int]
#    final_defect_type: DefectType
    final_vacancies: List[int]
    final_interstitials: List[int]
    defect_center_coord: Tuple[float, float, float]
    displacements: List[Displacement]
    neighboring_atom_indices: List[int]

    @property
    def symmetry_relation(self):
        return symmetry_relation(
            self.initial_point_group, self.final_point_group)

    @property
    def same_config_from_initial(self):
        return (self.initial_vacancies == self.final_vacancies
                and self.initial_interstitials == self.final_interstitials)

    def to(self, dir_name: Path = Path(".")):
        fd_local = Structure.from_sites(
            [self.final_structure[s] for s in self.neighboring_atom_indices])
        for v in self.final_vacancies:
            fd_local.append("X", coords=self.perfect_structure[v].frac_coords)
        fd_local.to(filename=dir_name / "POSCAR_final_local")

        neighbor_indices_in_p = [self.fd_to_p[i]
                                 for i in self.neighboring_atom_indices
                                 if self.fd_to_p[i] is not None]
        neighbor_indices_in_p.extend(self.final_vacancies)
        p_local = Structure.from_sites(
            [self.perfect_structure[s] for s in neighbor_indices_in_p])
        p_local.to(filename=dir_name / "POSCAR_perfect_local")

    def __repr__(self):
        lines = [f"Site symmetry: {self.initial_point_group} "
                 f"-> {self.final_point_group} ({self.symmetry_relation})",
                 f"Is same configuration: {self.same_config_from_initial}",
                 f"Drift distance: {self.drift_dist:5.3f}"]
        center_coords = [round(c, 3) for c in self.defect_center_coord]
        lines.append(f"Defect center: {center_coords}")

        lines.append("Removed atoms:")
        x = []
        for v in self.final_vacancies:
            site = self.perfect_structure[v]
            x.append([v, site.specie] + [round(fc, 2) for fc in site.frac_coords])
        lines.append(tabulate(x))

        lines.append("Added atoms:")
        x = []
        for i in self.final_interstitials:
            site = self.final_structure[i]
            x.append([i, site.specie] + [round(fc, 2) for fc in site.frac_coords])
        lines.append(tabulate(x))

        if self.same_config_from_initial is False:
            lines.append("Initially removed atoms:")
            x = []
            for v in self.initial_vacancies:
                site = self.perfect_structure[v]
                x.append([v, site.specie] + [round(fc, 2) for fc in site.frac_coords])
            lines.append(tabulate(x))

            lines.append("Initially added atoms:")
            x = []
            for i in self.initial_interstitials:
                site = self.initial_structure[i]
                x.append([i, site.specie] + [round(fc, 2) for fc in site.frac_coords])
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


def symmetrize_defect_structure(structure_symmetrizer: StructureSymmetrizer,
                                anchor_atom_idx: Optional[int] = None,
                                anchor_atom_coord: Optional[np.ndarray] = None
                                ) -> Structure:
    result = structure_symmetrizer.structure.copy()

    logger.info(f"The symmetry is {structure_symmetrizer.point_group}")

    origin_shift = structure_symmetrizer.spglib_sym_data["origin_shift"]
    inv_trans_mat = inv(structure_symmetrizer.spglib_sym_data["transformation_matrix"])
    coords = structure_symmetrizer.spglib_sym_data["std_positions"]
    new_coords = []
    for i in range(len(result)):
        new_coords.append(np.dot(inv_trans_mat, (coords[i] - origin_shift)))

    if anchor_atom_idx:
        offset = new_coords[anchor_atom_idx] - anchor_atom_coord
    else:
        offset = np.array([0.0]*3)
    new_coords = np.array(new_coords) - offset
    for i, coords in zip(result, new_coords):
        i.frac_coords = coords % 1
    return result

