# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict, Optional, Set

import numpy as np
from monty.json import MSONable
from pydefect.defaults import defaults
from pymatgen.core import Structure, Element


class Distances:
    def __init__(self,
                 structure: Structure,
                 center_coord: np.array,
                 dist_tol: float = None):
        self.structure = structure
        self.coord = center_coord
        self.dist_tol = dist_tol or defaults.dist_tol

    def distances(self, remove_self=True, specie=None) -> List[float]:
        result = []
        lattice = self.structure.lattice
        for site in self.structure:
            if specie and Element(specie) != site.specie:
                result.append(float("inf"))
                continue
            distance, _ = \
                lattice.get_distance_and_image(site.frac_coords, self.coord)
            if remove_self and distance < 1e-5:
                continue
            result.append(distance)

        return result

    def atom_idx_at_center(self, specie: str) -> Optional[int]:
        distances = self.distances(remove_self=False, specie=specie)
        sorted_dists = sorted(distances)
        if sorted_dists[0] > self.dist_tol:
            return None
        return np.argmin(distances)

    @property
    def shortest_distance(self) -> float:
        return min(self.distances())

    def coordination(self, include_on_site=False, cutoff_factor=None
                     ) -> "Coordination":
        cutoff_factor = cutoff_factor or defaults.cutoff_distance_factor
        cutoff = self.shortest_distance * cutoff_factor
        elements = [element.specie.name for element in self.structure]
        e_d = zip(elements, self.distances(remove_self=False))

        unsorted_distances = defaultdict(list)
        neighboring_atom_indices = []
        for i, (element, distance) in enumerate(e_d):
            if distance < cutoff:
                if include_on_site or \
                        (include_on_site is False and distance > 1e-5):
                    unsorted_distances[element].append(round(distance, 2))
                    neighboring_atom_indices.append(i)

        distance_dict = {}
        for element, distances in unsorted_distances.items():
            distance_dict[element] = sorted(distances)

        return Coordination(distance_dict, round(cutoff, 3),
                            neighboring_atom_indices)


@dataclass
class Coordination(MSONable):
    distance_dict: Dict[str, List]
    cutoff: float
    neighboring_atom_indices: List[int]
