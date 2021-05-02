# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import math
import warnings
from typing import List, Tuple

import numpy as np
from pydefect.analyzer.defect_structure_comparator import \
    DefectStructureComparator
from pydefect.analyzer.defect_structure_info import Displacement, \
    DefectStructureInfo, unique_point_group
from pydefect.defaults import defaults
from pymatgen import PeriodicSite, Structure, Lattice
from vise.util.structure_symmetrizer import StructureSymmetrizer
from vise.util.typing import GenCoords


def folded_coords(site: PeriodicSite,
                  center: GenCoords
                  ) -> GenCoords:
    _, image = site.distance_and_image_from_frac_coords(center)
    return tuple(site.frac_coords - image)


def fold_coords_in_structure(structure: Structure, center: GenCoords) -> None:
    for site in structure:
        _, image = site.distance_and_image_from_frac_coords(center)
        site.frac_coords -= image


def calc_drift(perfect: Structure,
               defect: Structure,
               center: List[float],
               d_to_p: List[int]) -> Tuple[int, float, GenCoords]:
    distances = []
    for site in defect:
        distances.append(site.distance_and_image_from_frac_coords(center)[0])
    anchor_atom_idx = int(np.argmax(distances))
    p_anchor_atom_idx = d_to_p[anchor_atom_idx]
    d_site = defect[anchor_atom_idx]
    p_coords = perfect[p_anchor_atom_idx].frac_coords
    distance, image = d_site.distance_and_image_from_frac_coords(p_coords)
    drift_vector = d_site.frac_coords - p_coords - image
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

            ini_dist = np.linalg.norm(initial_pos_vec)
            angle = calc_disp_angle(disp_dist, disp_vec, ini_dist,
                                    initial_pos_vec)

            result.append(Displacement(specie=p_elem,
                                       original_pos=initial_pos,
                                       final_pos=final_pos,
                                       distance_from_defect=ini_dist,
                                       disp_vector=tuple(disp_vec),
                                       displace_distance=disp_dist,
                                       angle=angle))
    return result


def calc_disp_angle(disp_dist, disp_vec, ini_dist, initial_pos_vec):
    inner_prod = sum(initial_pos_vec * disp_vec)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # ignore "RuntimeWarning: invalid value encountered in double_scalars"
        cos = round(inner_prod / (ini_dist * disp_dist), 10)
    result = float(round(180 * (1 - np.arccos(cos) / np.pi), 1))
    if math.isnan(result):
        result = None
    return result


def remove_dot(x):
    return "".join([s for s in x if s != "."])


def make_defect_structure_info(perfect: Structure,
                               initial: Structure,
                               final: Structure,
                               dist_tol: float,
                               symprec: float,
                               init_site_sym: str,
                               neighbor_cutoff_factor: float = None
                               ) -> "DefectStructureInfo":
    cutoff = neighbor_cutoff_factor or defaults.cutoff_distance_factor

    final = final.copy()
    symmetrizer = StructureSymmetrizer(final, symprec)
    final_site_sym = symmetrizer.point_group

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
        neighbor_atom_indices=comp_w_perf.neighboring_atom_indices(cutoff),
        neighbor_cutoff_factor=cutoff,
        drift_vector=tuple(drift_vec),
        drift_dist=drift_d,
        center=tuple(center),
        displacements=displacements)