# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from itertools import cycle
from typing import List, Optional, Tuple

import plotly.graph_objects as go
from adjustText import adjust_text
from labellines import labelLines
from matplotlib import pyplot as plt
from pydefect.analyzer.defect_energy import DefectEnergySummary
from pydefect.defaults import defaults
from vise.util.matplotlib import float_to_int_formatter


class DefectEnergiesMplSettings:
    def __init__(self,
                 colors: Optional[List[str]] = None,
                 line_width: float = 1.0,
                 vline_width: float = 0.5,
                 vline_color: str = "black",
                 vline_style: str = "-.",
                 vline_alpha: float = 0.7,
                 circle_size: int = 15,
                 tick_label_size: Optional[int] = 12,
                 title_font_size: Optional[int] = 15,
                 label_font_size: Optional[int] = 15,
                 defect_name_size: Optional[int] = 12,
                 charge_size: Optional[int] = 12):
        self.colors = cycle(colors) if colors else defaults.defect_energy_colors
        self.line_width = line_width
        self.circle_size = circle_size
        self.tick_label_size = tick_label_size
        self.title_font_size = title_font_size
        self.label_font_size = label_font_size
        self.defect_name_size = defect_name_size
        self.charge_size = charge_size

        self.vline = {"linewidth": vline_width,
                      "color":     vline_color,
                      "linestyle": vline_style,
                      "alpha": vline_alpha}


class DefectEnergyPlotter:
    def __init__(self,
                 defect_energy_summary: DefectEnergySummary,
                 chem_pot_label: "str",
                 allow_shallow: bool,
                 with_corrections: bool,
                 name_style: Optional[str],
                 x_range: Optional[Tuple[float, float]] = None,
                 y_range: Optional[Tuple[float, float]] = None,
                 vline_threshold: float = 0.02,
                 x_unit: Optional[str] = "eV",
                 y_unit: Optional[str] = "eV",
                 **plot_settings):

        self._title = defect_energy_summary.latexified_title
        self._supercell_vbm = defect_energy_summary.supercell_vbm
        self._supercell_cbm = defect_energy_summary.supercell_cbm
        self._x_range = x_range or (0, defect_energy_summary.cbm)
        defect_energy_summary.e_min = self._x_range[0]
        defect_energy_summary.e_max = self._x_range[1]
        charge_energies = defect_energy_summary.charge_energies(
            chem_pot_label, allow_shallow, with_corrections, self._x_range,
            name_style)
        self._cross_points = charge_energies.cross_point_dicts
        self._y_range = y_range or charge_energies.energy_range(space=0.2)
        self._vline_threshold = vline_threshold
        self._x_unit = x_unit
        self._y_unit = y_unit


class DefectEnergyPlotlyPlotter(DefectEnergyPlotter):
    def __init__(self, **kwargs):
        super().__init__(name_style="plotly", **kwargs)

    def create_figure(self):
        fig = go.Figure()
        fig.update_layout(
            title=f"Defect formation energy in {self._title}",
            xaxis_title=f"Fermi level ({self._x_unit})",
            yaxis_title=f"Energy ({self._y_unit})",
            title_font_size=30,
            font_size=24,
            width=900, height=700)

        for name, cp in self._cross_points.items():
            xs, ys = cp.t_all_sorted_points
            fig.add_trace(go.Scatter(x=xs, y=ys, name=name,
                                     text=cp.charge_list,
                                     hovertemplate=
                                     f'{name}<br>' +
                                     'Charges %{text}<br>' +
                                     'Energy: %{y:.2f}',
                                     line_width=3,
                                     marker_size=15))

        fig["layout"]["xaxis"]["range"] = [self._x_range[0] - 0.1,
                                           self._x_range[1] + 0.1]
        fig["layout"]["yaxis"]["range"] = self._y_range

        kwargs = dict(y=self._y_range + [None] + self._y_range,
                      line=dict(width=2, dash="dash"),
                      line_color="black")
        name = "unitcell"
        fig.add_trace(go.Scatter(
            x=[self._x_range[0]]*2 + [None] + [self._x_range[1]]*2,
            name=name, showlegend=True, **kwargs))

        x = [None] * 5
        if self._supercell_vbm > self._vline_threshold:
            x[0], x[1] = self._supercell_vbm, self._supercell_vbm
        if self._supercell_cbm < self._x_range[1] - self._vline_threshold:
            x[3], x[4] = self._supercell_cbm, self._supercell_cbm

        if len(set(x)) > 1:
            kwargs["name"] = "supercell"
            kwargs["line_color"] = "blue"
            kwargs["line"]["dash"] = "dot"
            fig.add_trace(go.Scatter(x=x, showlegend=True, **kwargs))

        return fig


class DefectEnergyMplPlotter(DefectEnergyPlotter):
    def __init__(self,
                 label_line: bool = True,
                 add_charges: bool = True,
                 **kwargs):
        super().__init__(name_style="mpl", **kwargs)
        self._mpl_defaults = \
            kwargs.get("mpl_defaults", DefectEnergiesMplSettings())
        self._label_line = label_line
        self._add_charges = add_charges
        self.plt = plt
        self._texts = []

    def construct_plot(self):
        self._add_energies()
        self._add_band_edges()
        self._set_x_range()
        self._set_y_range()
        self._set_labels()
        self._set_title()
        self._set_formatter()

        if self._add_charges:
            adjust_text(self._texts, force_points=(1.0, 2.5))

        if self._label_line:
            labelLines(plt.gca().get_lines(),
                       align=False,
                       fontsize=self._mpl_defaults.defect_name_size)
        else:
            ax = self.plt.gca()
            ax.legend(bbox_to_anchor=(1, 0.5), loc='center left')

        self.plt.tight_layout()

    def _add_energies(self):
        for name, cp in self._cross_points.items():
            color = next(self._mpl_defaults.colors)
            self.plt.plot(*cp.t_all_sorted_points, color=color,
                          linewidth=self._mpl_defaults.line_width,
                          label=name)
            if cp.t_inner_cross_points:
                self.plt.scatter(*cp.t_inner_cross_points, marker="o",
                                 color=color, s=self._mpl_defaults.circle_size)
            if self._add_charges:
                self._texts.extend(
                    [self.plt.text(x, y, charge, color=color,
                                   fontsize=self._mpl_defaults.charge_size)
                     for charge, (x, y) in cp.annotated_charge_positions.items()])

    def _set_x_range(self):
        self.plt.xlim(self._x_range)

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
        if self._supercell_vbm > self._vline_threshold:
            self.plt.axvline(x=self._supercell_vbm,
                             **self._mpl_defaults.vline)
            plt.text(self._supercell_vbm, self._y_range[1], 'supercell VBM',
                     size=8, ha='center', va='center', rotation='vertical',
                     backgroundcolor='white')
        if self._supercell_cbm < self._x_range[1] - self._vline_threshold:
            self.plt.axvline(x=self._supercell_cbm,
                             **self._mpl_defaults.vline)
            plt.text(self._supercell_cbm, self._y_range[1], 'supercell',
                     size=8, ha='center', va='center', rotation='vertical',
                     backgroundcolor='white')
