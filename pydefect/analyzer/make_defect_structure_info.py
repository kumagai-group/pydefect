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
from pymatgen import PeriodicSite, Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer
from vise.util.typing import GenCoords


def folded_coords(site: PeriodicSite,
                  center: GenCoords
                  ) -> GenCoords:
    _, image = site.distance_and_image_from_frac_coords(center)
    return tuple(site.frac_coords - image)


class MakeDefectStructureInfo:
    def __init__(self,
                 perfect: Structure,
                 initial: Structure,
                 final: Structure,
                 symprec: float,
                 dist_tol: float,
                 neighbor_cutoff_factor: float = None):

        self.cutoff = neighbor_cutoff_factor or defaults.cutoff_distance_factor
        self.symprec = symprec
        self.perfect, self.initial, self.final = perfect, initial, final
        self.dist_tol = dist_tol

        assert perfect.lattice == initial.lattice == final.lattice
        self.lattice = perfect.lattice

        self._orig_comp = DefectStructureComparator(final, perfect, dist_tol)
        self._orig_center = self._orig_comp.defect_center_coord
        self._calc_drift()

        self.shifted_final = final.copy()
        self.center = tuple(self._orig_center - self._drift_vector)
        for site in self.shifted_final:
            site.frac_coords -= np.array(self._drift_vector)

        self.comp_w_perf = DefectStructureComparator(
            self.shifted_final, perfect, dist_tol)
        self.comp_w_init = DefectStructureComparator(
            self.shifted_final, initial, dist_tol)

        self.defect_structure_info = DefectStructureInfo(
            shifted_final_structure=self.shifted_final,
            initial_site_sym=self.initial_site_sym,
            final_site_sym=self.final_site_sym,
            site_diff=self.comp_w_perf.make_site_diff(),
            site_diff_from_initial=self.comp_w_init.make_site_diff(),
            symprec=symprec,
            dist_tol=dist_tol,
            anchor_atom_idx=self._anchor_atom_idx,
            neighbor_atom_indices=self._neighbor_atom_indices,
            neighbor_cutoff_factor=self.cutoff,
            drift_vector=self._drift_vector,
            drift_dist=self._drift_distance,
            center=self.center,
            displacements=self.calc_displacements())

    @property
    def _neighbor_atom_indices(self):
        return self.comp_w_perf.neighboring_atom_indices(self.cutoff)

    @property
    def initial_site_sym(self):
        return self._unique_point_group(self.initial)

    @property
    def final_site_sym(self):
        return self._unique_point_group(self.final)

    def _unique_point_group(self, structure):
        symmetrizer = StructureSymmetrizer(structure, self.symprec)
        return unique_point_group(symmetrizer.point_group)

    def _calc_drift(self) -> None:
        distances = []
        for site in self.final:
            distances.append(
                site.distance_and_image_from_frac_coords(self._orig_center)[0])

        self._anchor_atom_idx = int(np.argmax(distances))
        p_anchor_atom_idx = self._orig_comp.d_to_p[self._anchor_atom_idx]
        d_site = self.final[self._anchor_atom_idx]
        p_coords = self.perfect[p_anchor_atom_idx].frac_coords
        self._drift_distance, image = \
            d_site.distance_and_image_from_frac_coords(p_coords)
        self._drift_vector = tuple(d_site.frac_coords - p_coords - image)

    def calc_displacements(self):
        result = []
        atom_mapping = self.comp_w_init.atom_mapping
        for d in range(len(self.shifted_final)):
            if d not in atom_mapping:
                result.append(None)
            else:
                p = atom_mapping[d]
                p_elem = str(self.initial[p].specie)
                d_elem = str(self.shifted_final[d].specie)
                if p_elem != d_elem:
                    result.append(None)
                    continue
                initial_pos = folded_coords(self.initial[p], self.center)
                final_pos = folded_coords(self.shifted_final[d], initial_pos)

                initial_pos_vec = self.lattice.get_cartesian_coords(
                    np.array(initial_pos) - self.center)
                initial_dist = np.linalg.norm(initial_pos_vec)

                disp_dist, t = self.lattice.get_distance_and_image(
                    self.shifted_final.frac_coords[d], self.initial.frac_coords[p])
                disp_vec = self.lattice.get_cartesian_coords(
                    self.shifted_final.frac_coords[d] - self.initial.frac_coords[p] - t)

                angle = self.calc_disp_angle(disp_dist, disp_vec, initial_dist,
                                             initial_pos_vec)

                result.append(Displacement(specie=p_elem,
                                           original_pos=initial_pos,
                                           final_pos=final_pos,
                                           distance_from_defect=initial_dist,
                                           disp_vector=tuple(disp_vec),
                                           displace_distance=disp_dist,
                                           angle=angle))
        return result

    @staticmethod
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

