# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from vise.analyzer.vesta.vesta_file import VestaFile
from pydefect.defaults import defaults
from pymatgen.core import Structure, DummySpecies, Element
from vise.util.typing import GenCoords


class MakeDefectVestaFile:
    def __init__(self,
                 defect_structure_info: DefectStructureInfo,
                 cutoff: float = defaults.show_structure_cutoff,
                 min_displace_w_arrows: float = 0.1,
                 arrow_factor: float = 3.0,
                 title: str = None):
        self.defect_structure = defect_structure_info.shifted_final_structure
        self.lattice = self.defect_structure.lattice

        self.removed_by_sub = defect_structure_info.site_diff.removed_by_sub
        self.inserted_by_sub = defect_structure_info.site_diff.inserted_by_sub
        self.removed = defect_structure_info.site_diff.removed
        self.inserted = defect_structure_info.site_diff.inserted
        self.removed_from_init = defect_structure_info.site_diff_from_initial.removed
        self.is_complex = defect_structure_info.site_diff.is_complex_defect

        self.center = defect_structure_info.center
        self.displacements = defect_structure_info.displacements

        self.cutoff = cutoff
        self.min_displace_w_arrows = min_displace_w_arrows
        self.arrow_factor = arrow_factor

        self._make_local_species_info()
        self._initial_local = self.make_initial_structure()
        self._final_local = self.make_final_structure()

        boundary = sum([[x - 0.5, x + 0.5]
                        for x in defect_structure_info.center], [])

        self.initial_vesta = \
            VestaFile(self._initial_local, title, boundary=boundary, show_label=False)
        self.final_vesta = \
            VestaFile(self._final_local, title, vectors=self._vectors,
                      vector_colors=self._vector_colors, boundary=boundary)

    def make_initial_structure(self):
        result = Structure(self.lattice, self._species, self._initial_coords)
        for _, specie, coord in self.removed_by_sub + self.removed:
            result.append(specie, coords=coord, properties={"name": f"{specie}"})
        fold_coords_in_structure(result, self.center)
        return result

    def make_final_structure(self):
        result = Structure(self.lattice, self._species, self._final_coords,
                           site_properties={"name": self._site_info})
        # inserted
        for _, specie, coords in self.inserted:
            result.append(specie, coords, properties={"name": f"{specie}_i"})
        # substituted
        for (_, i_specie, i_fcoord), (_, r_specie, r_fcoord) \
                in zip(self.inserted_by_sub, self.removed_by_sub):
            result.append(i_specie, i_fcoord,
                          properties={"name": f"{i_specie}_{r_specie}"})
        # removed. DummySpecies come later.
        for _, specie, coords in self.removed:
            result.append(DummySpecies(), coords,
                          properties={"name": f"Va_{specie}"})

        for _, specie, coords in self.removed_from_init:
            result.append(DummySpecies(), coords,
                          properties={"name": f"Initial_{specie}"})

        if self.is_complex:
            result.append(DummySpecies(), self.center,
                          properties={"name": "center"})

        fold_coords_in_structure(result, self.center)
        return result

    def _make_local_species_info(self):
        self._species, self._initial_coords, self._final_coords = [], [], []
        self._site_info, self._vectors, self._vector_colors = [], {}, []

        specie_idx = 1
        for specie, disp in zip(self.defect_structure, self.displacements):
            if disp and disp.distance_from_defect < self.cutoff:
                self._species.append(Element(disp.specie))
                self._initial_coords.append(disp.original_pos)
                self._final_coords.append(specie.frac_coords)
                info = f"{round(disp.distance_from_defect, 1)}"

                if disp and disp.displace_distance > self.min_displace_w_arrows:
                    self._vectors[specie_idx] = \
                        (v * self.arrow_factor for v in disp.disp_vector)
                    # max=200 instead of 255 since (255, 255, 255) is white.
                    rgb = int(200 * disp.angle / 180)
                    self._vector_colors.append((rgb, rgb, 100))
                    info += f"_{round(disp.displace_distance, 2)}"

                self._site_info.append(info)
                specie_idx += 1


def fold_coords_in_structure(structure: Structure, center: GenCoords) -> None:
    for site in structure:
        _, image = site.distance_and_image_from_frac_coords(center)
        site.frac_coords -= image