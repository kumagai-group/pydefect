# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import cycle

from monty.design_patterns import singleton
from vise.defaults import DefaultsBase


@singleton
class Defaults(DefaultsBase):
    def __init__(self):
        self._symmetry_length_tolerance = 0.1
        self._symmetry_angle_tolerance = 5.0
        self._ewald_accuracy = 15.0
        self._e_above_hull = 1e-5
        self._cutoff_distance_factor = 1.3
        self._show_structure_cutoff = 5.0
        self._displace_distance = 0.2
        self._dist_tol = 1.0
        self._ele_neg_diff = 2.0
        self._similar_orb_criterion = 0.2
        self._similar_energy_criterion = 0.5  # in eV
        self._state_occupied_threshold = 0.07
        self._eigval_range = 1.0
        self._defect_energy_colors = \
            ["xkcd:blue", "xkcd:brown", "xkcd:crimson", "xkcd:darkgreen",
             "xkcd:gold", "xkcd:magenta", "xkcd:orange", "xkcd:darkblue",
             "xkcd:navy", "xkcd:red", "xkcd:olive", "xkcd:black", "xkcd:indigo"]
        self._abs_strange_energy = 100.0
        self._localized_orbital_radius = 3.0
        self._localized_orbital_fraction_wrt_uniform = 0.7

        self.set_user_settings(yaml_filename="pydefect.yaml")

    @property
    def symmetry_length_tolerance(self):
        return self._symmetry_length_tolerance

    @property
    def symmetry_angle_tolerance(self):
        return self._symmetry_angle_tolerance

    @property
    def ewald_accuracy(self):
        return self._ewald_accuracy

    @property
    def e_above_hull(self):
        return self._e_above_hull

    @property
    def cutoff_distance_factor(self):
        return self._cutoff_distance_factor

    @property
    def show_structure_cutoff(self):
        return self._show_structure_cutoff

    @property
    def displace_distance(self):
        return self._displace_distance

    @property
    def dist_tol(self):
        return self._dist_tol

    @property
    def ele_neg_diff(self):
        return self._ele_neg_diff

    @property
    def similar_orb_criterion(self):
        return self._similar_orb_criterion

    @property
    def similar_energy_criterion(self):
        return self._similar_energy_criterion

    @property
    def state_occupied_threshold(self):
        return self._state_occupied_threshold

    @property
    def state_unoccupied_threshold(self):
        return 1 - self._state_occupied_threshold

    @property
    def eigval_range(self):
        return self._eigval_range

    @property
    def defect_energy_colors(self):
        return cycle(self._defect_energy_colors)

    @property
    def abs_strange_energy(self):
        return self._abs_strange_energy

    @property
    def localized_orbital_radius(self):
        return self._localized_orbital_radius

    @property
    def localized_orbital_fraction_wrt_uniform(self):
        return self._localized_orbital_fraction_wrt_uniform


defaults = Defaults()
