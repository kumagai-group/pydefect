# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import cycle
from typing import List, Optional

from matplotlib import pyplot as plt
from vise.util.matplotlib import float_to_int_formatter

from pydefect.analyzer.defect_energy import DefectEnergy
from pydefect.defaults import defaults


class DefectEnergiesMplSettings:
    def __init__(self,
                 colors: Optional[List[str]] = None,
                 line_width: float = 1.0,
                 band_edge_line_width: float = 1.0,
                 band_edge_line_color: str = "black",
                 band_edge_line_style: str = "-",
                 supercell_band_edge_line_width: float = 0.5,
                 supercell_band_edge_line_color: str = "black",
                 supercell_band_edge_line_style: str = "-.",
                 circle_size: int = 15,
                 tick_label_size: Optional[int] = 12,
                 title_font_size: Optional[int] = 15,
                 label_font_size: Optional[int] = 15):
        self.colors = cycle(colors) if colors else defaults.defect_energy_colors
        self.line_width = line_width
        self.circle_size = circle_size
        self.tick_label_size = tick_label_size
        self.title_font_size = title_font_size
        self.label_font_size = label_font_size

        self.vline = {"linewidth": band_edge_line_width,
                      "color": band_edge_line_color,
                      "linestyle": band_edge_line_style}

        self.svline = {"linewidth": supercell_band_edge_line_width,
                       "color": supercell_band_edge_line_color,
                       "linestyle": supercell_band_edge_line_style}


class DefectEnergyPlotter:
    def __init__(self,
                 title: str,
                 defect_energies: List[DefectEnergy],
                 vbm: float,
                 cbm: float,
                 supercell_vbm: float,
                 supercell_cbm: float,
                 y_range: Optional[List[float]] = None,
                 x_unit: Optional[str] = "eV",
                 y_unit: Optional[str] = "eV",
                 mpl_defaults: Optional[DefectEnergiesMplSettings]
                 = DefectEnergiesMplSettings()
                 ):

        self._title = title
        self._defect_energies = defect_energies
        self._vbm = vbm
        self._cbm = cbm
        self._supercell_vbm = supercell_vbm
        self._supercell_cbm = supercell_cbm
        self._y_range = y_range

        self._x_unit = x_unit
        self._y_unit = y_unit

        self._mpl_defaults = mpl_defaults
        self.plt = plt

    def construct_plot(self):
        self._add_energies()
        self._add_band_edges()
        self._set_x_range()
        self._set_y_range()
        self._set_labels()
        self._set_title()
        self._set_formatter()
        self.plt.tight_layout()
        self.plt.gca().legend(bbox_to_anchor=(1, 0.5), loc='center left')

    def _add_energies(self):
        for de in self._defect_energies:
            color = next(self._mpl_defaults.colors)
            cp = de.cross_points(self._supercell_vbm, self._supercell_cbm)
            self.plt.plot(*cp.t_all_sorted_points, color=color,
                          linewidth=self._mpl_defaults.line_width,
                          label=de.name)
            if cp.t_inner_cross_points:
                self.plt.scatter(*cp.t_inner_cross_points, marker="o",
                                 color=color, s=self._mpl_defaults.circle_size)

    def _set_x_range(self):
        self.plt.xlim(self._supercell_vbm, self._supercell_cbm)

    def _set_y_range(self):
        if self._y_range:
            self.plt.ylim(self._y_range[0], self._y_range[1])

    def _set_labels(self):
        self.plt.xlabel(f"Fermi level ({self._x_unit})",
                        size=self._mpl_defaults.label_font_size)
        self.plt.ylabel(f"Energy ({self._y_unit})",
                        size=self._mpl_defaults.label_font_size)

    def _set_title(self):
        self.plt.title(self._title, size=self._mpl_defaults.title_font_size)

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.yaxis.set_major_formatter(float_to_int_formatter)
        axis.tick_params(labelsize=self._mpl_defaults.tick_label_size)

    def _add_band_edges(self):
        self.plt.axvline(x=self._vbm, **self._mpl_defaults.vline)
        self.plt.axvline(x=self._cbm, **self._mpl_defaults.vline)
        self.plt.axvline(x=self._supercell_vbm, **self._mpl_defaults.svline)
        self.plt.axvline(x=self._supercell_cbm, **self._mpl_defaults.svline)

