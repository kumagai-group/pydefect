# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from itertools import product
from typing import Tuple, List

import numpy as np
from pydefect.cli.vasp.make_efnv_correction import calc_max_sphere_radius
from pymatgen import Lattice
from pymatgen.io.vasp import Chgcar


@dataclass
class Grids:
    lattice: Lattice
    dim: Tuple[int, int, int]
    distance_data: np.ndarray

    def dump(self, filename="grids.npz"):
        np.savez(filename, matrix=self.lattice.matrix,
                 distance_data=self.distance_data)

    @classmethod
    def from_file(cls, filename="grids.npz"):
        loaded_dict = np.load(filename)
        lattice = Lattice(loaded_dict["matrix"])
        return cls(dim=loaded_dict["distance_data"].shape, lattice=lattice,
                   distance_data=loaded_dict["distance_data"])

    @classmethod
    def from_chgcar(cls, chgcar: Chgcar):
        lattice, dim = chgcar.structure.lattice, chgcar.dim
        grid_points = [[x / dim[0], y / dim[1], z / dim[2]]
                       for (x, y, z) in product(*[list(range(i)) for i in dim])]
        # Use boolean indexing to find charges within the desired distance.
        # data[:, 0]: shifted_coords
        # data[:, 1]: distances
        # data[:, 2]: sequential indices
        # data[:, 3]: images
        distances_data = np.array(lattice.get_all_distances(
            grid_points, (0.0, 0.0, 0.0)))
        return cls(lattice, dim, distances_data[:, 0].reshape(dim))

    def shifted_distance_data(self, center: List[int]):
        return np.roll(np.roll(np.roll(self.distance_data, center[0], axis=0),
                               center[1], axis=1), center[2], axis=2)

    def spherical_dist(self,
                       data: np.ndarray,
                       center: List[int],
                       distance_bins: np.ndarray):
        assert distance_bins[-1] <= calc_max_sphere_radius(self.lattice.matrix)
        shifted_dist = self.shifted_distance_data(center)
        _sum, _ = np.histogram(shifted_dist, distance_bins, weights=data)
        counts, _ = np.histogram(shifted_dist, distance_bins)
        histogram = _sum / counts / self.lattice.volume
        return histogram.tolist()

#        num_bins = int(np.ceil(radius / bin_interval))
#        distance_bins = [bin_interval * i for i in range(num_bins)] + [radius]
    # distance_bins: List[float]
