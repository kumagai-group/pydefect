# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import groupby, cycle
from typing import Optional, List

import plotly.graph_objects as go
from matplotlib import pyplot as plt
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.corrections.gkfo_correction import GkfoCorrection
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
                 sites: List[PotentialSite],
                 defect_region_radius: float,
                 ave_pot_diff: float,
                 x_unit: Optional[str] = "Å",
                 y_unit: Optional[str] = "V"):
        self._title = title
        self.sites = sites
        self.defect_region_radius = defect_region_radius
        self.ave_pot_diff = ave_pot_diff

        self._x_unit = x_unit
        self._y_unit = y_unit

    @classmethod
    def from_efnv_corr(cls,
                       title,
                       efnv_correction: ExtendedFnvCorrection,
                       x_unit: Optional[str] = "Å",
                       y_unit: Optional[str] = "V"):
        return cls(title=title,
                   sites=efnv_correction.sites,
                   defect_region_radius=efnv_correction.defect_region_radius,
                   ave_pot_diff=efnv_correction.average_potential_diff,
                   x_unit=x_unit, y_unit=y_unit)

    @classmethod
    def from_gkfo_corr(cls,
                       title,
                       gkfo_correction: GkfoCorrection,
                       x_unit: Optional[str] = "Å",
                       y_unit: Optional[str] = "V"):
        efnv = gkfo_correction.init_efnv_correction
        ave_pot_diff = gkfo_correction.average_potential_diff_by_addition
        return cls(title=title,
                   sites=gkfo_correction.gkfo_sites,
                   defect_region_radius=efnv.defect_region_radius,
                   ave_pot_diff=ave_pot_diff,
                   x_unit=x_unit, y_unit=y_unit)


class SitePotentialPlotlyPlotter(SitePotentialPlotter):
    def create_figure(self):
        common = {"mode": "markers", "marker_size": 10}
        fig = go.Figure()
        fig.update_layout(
            title="Potential",
            xaxis_title=f"Distance ({self._x_unit})",
            yaxis_title=f"Potential ({self._y_unit})",
            title_font_size=30,
            font_size=24, width=800, height=700)
        all_distances = []
        pc_pot_distances = []
        all_potentials = []
        pc_potentials = []
        potential_diffs = []
        key = lambda x: x.specie
        for _, grouped_sites in groupby(sorted(self.sites, key=key), key=key):
            distances = []
            potentials = []
            atom_idxs = []
            for site in grouped_sites:
                distances.append(site.distance)
                all_distances.append(site.distance)
                all_potentials.append(site.potential)
                potentials.append(site.potential)
#                atom_idxs.append(site.)

                if site.pc_potential:
                    pc_pot_distances.append(site.distance)
                    pc_potentials.append(site.pc_potential)
                    potential_diffs.append(site.diff_pot)

            fig.add_trace(go.Scatter(x=distances, y=potentials,
                                     text=atom_idxs,
                                     name=site.specie,
                                     marker_symbol="circle",
                                     **common))

        min_y = min(all_potentials + pc_potentials + potential_diffs) - 0.4
        max_y = max(all_potentials + pc_potentials + potential_diffs) + 0.4
        max_x = max(all_distances) + 0.4
        fig.add_shape(dict(type="line",
                           x0=self.defect_region_radius,
                           x1=self.defect_region_radius,
                           y0=min_y, y1=max_y,
                           line=dict(color="royalblue", width=1.5, dash="dot")))
        fig["layout"]["xaxis"]["range"] = [0, max_x]
        fig["layout"]["yaxis"]["range"] = [min_y, max_y]

        fig.add_trace(go.Scatter(x=pc_pot_distances, y=pc_potentials,
                                 name="PC", marker_symbol="diamond",
                                 **common))
        fig.add_trace(go.Scatter(x=pc_pot_distances, y=potential_diffs,
                                 name="Diff", marker_symbol="cross",
                                 **common))

        fig.add_shape(dict(type="line",
                           x0=self.defect_region_radius, x1=max_x,
                           y0=self.ave_pot_diff, y1=self.ave_pot_diff,
                           line=dict(color="red", width=1.5, dash="dot"),
                           name="alignment"))

        return fig


class SitePotentialMplPlotter(SitePotentialPlotter):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._mpl_defaults = kwargs.get("mpl_defaults",
                                        PotentialPlotterMplSettings())
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
        pc_pot_distances = []
        pc_potentials = []
        potential_diffs = []
        key = lambda x: x.specie
        for _, grouped_sites in groupby(sorted(self.sites, key=key), key=key):
            distances = []
            potentials = []
            for site in grouped_sites:
                distances.append(site.distance)
                all_distances.append(site.distance)
                potentials.append(site.potential)

                if site.pc_potential:
                    pc_pot_distances.append(site.distance)
                    pc_potentials.append(site.pc_potential)
                    potential_diffs.append(site.diff_pot)

            color = next(self._mpl_defaults.colors)
            self.plt.scatter(distances, potentials, marker="o", color=color,
                             label=str(site.specie))

        self._max_distance = max(all_distances)

        self.plt.scatter(pc_pot_distances, pc_potentials, marker="1",
                         color="b", label="point charge")
        self.plt.scatter(pc_pot_distances, potential_diffs, marker="+",
                         color="r", label="potential difference")

        self.plt.legend(bbox_to_anchor=(0, 1.05), loc='upper left',
                        borderaxespad=0, fontsize=8)

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