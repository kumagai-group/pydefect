# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import cycle
from pathlib import Path

from monty.design_patterns import singleton
from vise.user_settings import UserSettings


@singleton
class Defaults:
    def __init__(self):
        self._symmetry_length_tolerance = 0.1
        self._symmetry_angle_tolerance = 5.0
        self._ewald_accuracy = 15.0
        self._e_above_hull = 1e-5
        self._e_above_hull = 1e-5
        self._cutoff_distance_factor = 1.3
        self._displace_distance = 0.2
        self._same_distance_criterion = 0.1
        self._kpoint_density = 5.0
        self._band_gap_criterion = 0.2  # in eV
        self._similar_orb_criterion = 0.12
        self._similar_energy_criterion = 0.2  # in eV
        self._eigval_range = 1.0
        self._integer_criterion = 0.1
        self._defect_energy_colors = \
            ["xkcd:blue", "xkcd:brown", "xkcd:crimson", "xkcd:darkgreen",
             "xkcd:gold", "xkcd:magenta", "xkcd:orange", "xkcd:darkblue",
             "xkcd:navy", "xkcd:red", "xkcd:olive", "xkcd:black", "xkcd:indigo"]
        self._outcar = "OUTCAR"
        self._contcar = "CONTCAR"
        self._vasprun = "vasprun.xml"
        self._procar = "PROCAR"

        user_settings = UserSettings(yaml_filename="pydefect.yaml")
        self.yaml_files = user_settings.yaml_files_from_root_dir
        self.user_settings = user_settings.user_settings

        for k, v in self.user_settings.items():
            if hasattr(self, k):
                self.__setattr__("_" + k, v)

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
    def displace_distance(self):
        return self._displace_distance

    @property
    def same_distance_criterion(self):
        return self._same_distance_criterion

    @property
    def kpoint_density(self):
        return self._kpoint_density

    @property
    def similar_orb_criterion(self):
        return self._similar_orb_criterion

    @property
    def similar_energy_criterion(self):
        return self._similar_energy_criterion

    @property
    def eigval_range(self):
        return self._eigval_range

    @property
    def band_gap_criterion(self):
        return self._band_gap_criterion

    @property
    def integer_criterion(self):
        return self._integer_criterion

    @property
    def defect_energy_colors(self):
        return cycle(self._defect_energy_colors)

    @property
    def outcar(self):
        return Path(self._outcar)

    @property
    def contcar(self):
        return Path(self._contcar)

    @property
    def vasprun(self):
        return Path(self._vasprun)

    @property
    def procar(self):
        return Path(self._procar)


defaults = Defaults()
