# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from itertools import product
from typing import List

import numpy as np
from pydefect.cli.vasp.make_efnv_correction import calc_max_sphere_radius
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Chgcar


class RadialDist:
    def __init__(self, parchg: Chgcar, center_coords: List[float]):
        self.data = parchg.spin_data
        self.dim = parchg.dim
        self.lattice = parchg.structure.lattice
        self.center = center_coords
        self.radius = calc_max_sphere_radius(self.lattice.matrix)
        self._distances_data = None  # lazy evaluation

    @property
    def distances_data(self):
        if self._distances_data is None:
            grid_points = [[x / self.dim[0], y / self.dim[1], z / self.dim[2]]
                           for (x, y, z) in
                           product(*[list(range(i)) for i in self.dim])]

            # Use boolean indexing to find charges within the desired distance.
            # data[:, 0]: shifted_coords
            # data[:, 1]: distances
            # data[:, 2]: sequential indices
            # data[:, 3]: images
            self._distances_data = np.array(self.lattice.get_points_in_sphere(
                grid_points, self.center, self.radius))

        return self._distances_data

    def histogram(self, spin: Spin, nbins: int = 15):
        inds = self.distances_data[:, 1] <= self.radius
        dists = self.distances_data[inds, 1]
        data_inds = np.rint(np.mod(list(self.distances_data[inds, 0]), 1) *
                            np.tile(self.dim, (len(dists), 1))).astype(int)

        vals = [self.data[spin][x, y, z] for x, y, z in data_inds]

        hist, edges = np.histogram(dists, nbins,
                                   range=[0, self.radius], weights=vals)
        hist_numbers, _ = np.histogram(dists, nbins, range=[0, self.radius])
        mesh_distance = edges[1] - edges[0]

        hist_data = np.zeros((nbins, 2))
        hist_data[:, 0] = [sum(edges[i:i + 2]) / 2 for i in range(nbins)]

        # 4pi * r^2 * rho
        integrated_volume = 4 * np.pi * hist_data[:, 0] ** 2
        density = hist / hist_numbers / self.lattice.volume
        hist_data[:, 1] = integrated_volume * density

        half_point = None
        for i in range(hist_data[:, 1].size):
            if sum(hist_data[:i + 1, 1]) > 0.5:
                # Obtain from the calculation of the area of a trapezoid
                x0 = hist_data[i - 1, 0]
                y0 = hist_data[i, 1]
                y1 = hist_data[i + 1, 1]
                half_point = (0.5 - sum(hist_data[:i, 1]) * mesh_distance) \
                    * 2 / (y0 + y1) + x0
                break

        summed = sum(hist_data[:, 1]) * mesh_distance
        return hist_data, half_point, summed


