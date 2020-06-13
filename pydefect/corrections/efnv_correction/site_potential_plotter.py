# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import groupby, cycle
from typing import Optional, List

from matplotlib import pyplot as plt
from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection
from pydefect.defaults import defaults
from vise.util.matplotlib import float_to_int_formatter


class PotentialPlotterMplSettings:
    def __init__(self,
                 colors: Optional[List[str]] = None,
                 line_width: float = 1.0,
                 radius_line_width: float = 1.0,
                 radius_line_color: str = "black",
                 radius_line_style: str = "-.",
                 ave_pot_diff_line_width: float = 1.5,
                 ave_pot_diff_line_color: str = "red",
                 ave_pot_diff_line_style: str = ":",
                 circle_size: int = 15,
                 zero_line_width: float = 0.5,
                 zero_line_color: str = "black",
                 zero_line_style: str = ":",
                 tick_label_size: Optional[int] = 12,
                 title_font_size: Optional[int] = 15,
                 label_font_size: Optional[int] = 15):
        self.colors = cycle(colors) if colors else defaults.defect_energy_colors
        self.line_width = line_width
        self.circle_size = circle_size
        self.tick_label_size = tick_label_size
        self.title_font_size = title_font_size
        self.label_font_size = label_font_size

        self.vline = {"linewidth": radius_line_width,
                      "color": radius_line_color,
                      "linestyle": radius_line_style}

        self.hline = {"linewidth": ave_pot_diff_line_width,
                      "color": ave_pot_diff_line_color,
                      "linestyle": ave_pot_diff_line_style}

        self.zero_line = {"linewidth": zero_line_width,
                          "color": zero_line_color,
                          "linestyle": zero_line_style}


class SitePotentialPlotter:
    def __init__(self,
                 title: str,
                 efnv_correction: ExtendedFnvCorrection,
                 x_unit: Optional[str] = "Å",
                 y_unit: Optional[str] = "V",
                 mpl_defaults: Optional[PotentialPlotterMplSettings]
                 = PotentialPlotterMplSettings()):

        self._title = title
        self.sites = efnv_correction.sites
        self.defect_region_radius = efnv_correction.defect_region_radius
        self.ave_pot_diff = efnv_correction.average_potential_diff

        self._x_unit = x_unit
        self._y_unit = y_unit

        self._mpl_defaults = mpl_defaults
        self.plt = plt

    def construct_plot(self):
        self._add_potentials()
        self._add_lines()
        self._set_labels()
        self._set_x_range()
        self._set_title()
        self._set_formatter()
        self.plt.tight_layout()
#        self.plt.gca().legend(bbox_to_anchor=(1, 0.5), loc='center left')

    def _add_potentials(self):
        all_distances = []
        pc_potentials = []
        potential_diffs = []
        for _, grouped_sites in groupby(self.sites, key=lambda x: x.specie):
            distances = []
            potentials = []
            for site in grouped_sites:
                distances.append(site.distance)
                potentials.append(site.potential)

                all_distances.append(site.distance)
                pc_potentials.append(site.pc_potential)
                potential_diffs.append(site.diff_pot)

            color = next(self._mpl_defaults.colors)
            self.plt.scatter(distances, potentials, marker="o", color=color)

        self._max_distance = max(all_distances)

        self.plt.scatter(all_distances, pc_potentials, marker="1", color="b")
        self.plt.scatter(all_distances, potential_diffs, marker="+", color="r")

    def _add_lines(self):
        self.plt.axvline(x=self.defect_region_radius,
                         **self._mpl_defaults.vline)
        self.plt.plot([self.defect_region_radius, self._max_distance * 1.1],
                      [self.ave_pot_diff, self.ave_pot_diff],
                      **self._mpl_defaults.hline)
        self.plt.axhline(y=0, **self._mpl_defaults.zero_line)

    def _set_x_range(self):
        self.plt.xlim(0, self._max_distance * 1.05)

    def _set_title(self):
        self.plt.title(self._title, size=self._mpl_defaults.title_font_size)

    def _set_labels(self):
        self.plt.xlabel(f"Distance from a defect ({self._x_unit})",
                        size=self._mpl_defaults.label_font_size)
        self.plt.ylabel(f"Potential ({self._y_unit})",
                        size=self._mpl_defaults.label_font_size)

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.xaxis.set_major_formatter(float_to_int_formatter)
        axis.yaxis.set_major_formatter(float_to_int_formatter)
        axis.tick_params(labelsize=self._mpl_defaults.tick_label_size)