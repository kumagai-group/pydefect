# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import math
from dataclasses import dataclass
from typing import Optional, List, Tuple

import numpy as np
from monty.json import MSONable
from numpy.linalg import inv
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Structure, Lattice
from vise.util.enum import ExtendedEnum
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


def fold_coords(structure: Structure, center: List[float]):
    for site in structure:
        _, image = site.distance_and_image_from_frac_coords(center)
        site.frac_coords -= image


def calc_drift_dist(perfect: Structure, defect: Structure, center: List[float],
                    d_to_p: List[int]):
    distances = []
    for site in defect:
        distances.append(site.distance_and_image_from_frac_coords(center)[0])
    anchor_atom_idx = np.argmax(distances)
    p_anchor_atom_idx = d_to_p[anchor_atom_idx]
    return defect[anchor_atom_idx].distance(perfect[p_anchor_atom_idx])


def calc_displacements(perfect: Structure, defect: Structure,
                       center: List[float], d_to_p: List[int]):
    result = []
    lattice: Lattice = defect.lattice
    for d, p in enumerate(d_to_p):
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

            if disp_dist > 0.1 and 0 <= angle < 15:
                annotation = "inward"
            elif disp_dist > 0.1 and 75 <= angle < 90:
                annotation = "outward"
            else:
                annotation = None
            result.append(Displacement(specie=elem,
                                       original_pos=orig_pos,
                                       final_pos=final_pos,
                                       distance_from_defect=round(ini_dist, 3),
                                       displace_distance=round(disp_dist, 3),
                                       angle=angle,
                                       annotation=annotation))
    return result


def make_defect_structure_info(perfect: Structure,
                               defect_entry: DefectEntry,
                               calc_results: CalcResults):
    initial = defect_entry.structure.copy()
    final = calc_results.structure.copy()
    assert perfect.lattice == initial.lattice == final.lattice

    comp_i = DefectStructureComparator(initial, perfect)
    comp_f = DefectStructureComparator(final, perfect)

    center = comp_f.defect_center_coord

    fold_coords(initial, center)
    fold_coords(final, center)

    drift_dist = calc_drift_dist(perfect, final, center, comp_f.d_to_p)
    assert drift_dist < 0.01

    displacements = calc_displacements(perfect, final, center, comp_f.d_to_p)

    return DefectStructureInfo(initial_point_group=defect_entry.site_symmetry,
                               final_point_group=calc_results.site_symmetry,
                               initial_structure=initial,
                               final_structure=final,
                               perfect_structure=perfect,
                               drift_dist=drift_dist,
                               initial_vacancies=comp_i.vacant_indices,
                               initial_interstitials=comp_i.inserted_indices,
                               final_vacancies=comp_f.vacant_indices,
                               final_interstitials=comp_f.inserted_indices,
                               defect_center_coord=tuple(center),
                               displacements=displacements)


class DefectType(MSONable, ExtendedEnum):
    vacancy = "vacancy"
    interstitial = "interstitial"
    substituted = "substituted"
    complex = "complex"
#    split = "split"


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
    annotation: str = None

    # def __repr__(self):
    #     return f" {self.specie:>2} {self.site:>3} {self.initial_distance:.2f>5} " \
    #            f"{self.displace_distance:} {} -> {}"


@dataclass
class DefectStructureInfo(MSONable):
    initial_point_group: str
    final_point_group: str
    initial_structure: Structure  # frac coords are shifted.
    final_structure: Structure  # frac coords are shifted.
    perfect_structure: Structure
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
#    neighboring_atom_indices: List[int]

    def symmetry_relation(self):
        return
    # def __repr__(self):
    #     lines = [f"Site symmetry: {self.initial_sym} -> {self.final_sym}",
    #              f"Transferred atoms:"]
    #     for v in self.comparator.vacant_indices:
    #         site = self.perfect_structure[v]
    #         coords = [round(fc, 2) for fc in site.frac_coords]
    #         lines.append(f" - {v:>3} {site.specie:>2} {coords}")
    #     for i in self.comparator.inserted_indices:
    #         site = self.defect_structure[i]
    #         coords = [round(fc, 2) for fc in site.frac_coords]
    #         lines.append(f" + {i:>3} {site.specie:>2} {coords}")

        # lines.append("Displacements:")
        # return "\n".join(lines)


def fold_frac_coords(frac_coords, center=None):
    center = center if center is not None else [0., 0., 0.]
    for c in center:
        assert 0.0 <= c < 1.0
    result = []
    for f, c in zip(frac_coords, center):
        f = f % 1
        if f < c:
            result.append(f if c - f <= f + 1.0 - c else f + 1.0)
        else:
            result.append(f if f - c <= c - (f - 1.0) else f - 1.0)
    return np.array(result)


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

