# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from collections import defaultdict
from typing import List, Dict

import numpy as np
from pymatgen import Structure


class Distances:
    def __init__(self,
                 structure: Structure,
                 coord: np.array):
        self.structure = structure
        self.coord = coord

    def distances(self, remove_self=True) -> List[float]:
        result = []
        for frac_coord in self.structure.frac_coords:
            lattice = self.structure.lattice
            distance, _ = lattice.get_distance_and_image(frac_coord, self.coord)
            if remove_self and distance < 1e-5:
                continue
            result.append(distance)

        return result

    @property
    def shortest_distance(self) -> float:
        return min(self.distances())

    def coord_distances(self, cutoff_distance_factor: float
                        ) -> Dict[str, List[float]]:
        cutoff = self.shortest_distance * cutoff_distance_factor

        elements = [element.specie.name for element in self.structure]
        e_d = zip(elements, self.distances(False))

        unsorted_distances = defaultdict(list)
        for element, distance in filter(lambda x: 1e-5 < x[1] < cutoff, e_d):
            unsorted_distances[element].append(round(distance, 2))

        result = {}
        for element, distances in unsorted_distances.items():
            result[element] = sorted(distances)

        return result
