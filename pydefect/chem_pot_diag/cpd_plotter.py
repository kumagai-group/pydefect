# -*- coding: utf-8 -*-
from abc import abstractmethod, ABC
from itertools import cycle
from typing import Optional, List

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag
from vise.util.plotly_util import sort_coords, make_triangles
from vise.util.string import latexify


class CpdMplSettings:
    def __init__(self):
        # Avoid type 3 font, http://phyletica.org/matplotlib-fonts/
        # https://stackoverflow.com/questions/34387893/output-matplotlib-figure-to-svg-with-text-as-text-not-curves
        self.rc_params = {'pdf.fonttype':  42, 'ps.fonttype': 42,
                          'svg.fonttype': "none"}

        self.grid_2d = {"color": 'black',
                        "alpha": 0.8,  # grid opacity
                        "linestyle": 'dashed',
                        "linewidth": 0.5}

        # Default zorder
        # Patch / PatchCollection	1
        # Line2D / LineCollection	2
        # Text	3
        # Big number of zorder is needed as it increases with num of simplexes
        self.label = {"size": 'smaller', "zorder": 100}

        self.alphabet_label = {"color": "blue", "weight": "bold", "zorder": 101}

        self.target_comp_text_color = "blue"
        self.comp_text_color = "black"

        self.simplex_2d = {"color": "black", "linewidth": 2, "zorder": 1}


class ChemPotDiagMplPlotter(ABC):
    def __init__(self,
                 cpd: ChemPotDiag,
                 mpl_defaults: Optional[CpdMplSettings] = CpdMplSettings()):
        self.cpd = cpd
        self._mpl_defaults = mpl_defaults
        self._add_ax()

    @abstractmethod
    def _add_ax(self):
        pass

    def draw_diagram(self, title: Optional[str] = None) -> plt:
        plt.rcParams.update(self._mpl_defaults.rc_params)
        self._plot()
        self._add_alphabet_labels()
        self._set_title(title)
        self._set_lims()
        self._set_labels()
        self._set_grid()
        plt.tight_layout()

        return plt

    @abstractmethod
    def _set_lims(self):
        pass

    @abstractmethod
    def _set_labels(self):
        pass

    @abstractmethod
    def _set_grid(self):
        pass

    def _plot(self):
        for comp, position in self.cpd.comp_centers.items():
            self.draw_simplex(comp)
            self._ax.text(*position,
                          latexify(comp),
                          color=self._text_color(comp),
                          **self._mpl_defaults.label,
                          **self._text_kwargs)

    @abstractmethod
    def draw_simplex(self, formula):
        pass

    @property
    @abstractmethod
    def _text_kwargs(self):
        pass

    def _text_color(self, formula):
        if self.cpd.target_coords and formula == self.cpd.target:
            return self._mpl_defaults.target_comp_text_color
        return self._mpl_defaults.comp_text_color

    def _set_title(self, title):
        if not title:
            title = f"Chemical potential diagram of {self.cpd.chemical_system}"
        self._ax.set_title(title)

    def _add_alphabet_labels(self):
        if self.cpd.target:
            for label, vertex_coords in self.cpd.target_coords.items():
                self._ax.text(*transpose(vertex_coords), label,
                              **self._mpl_defaults.alphabet_label)


class ChemPotDiag2DMplPlotter(ChemPotDiagMplPlotter):

    def _add_ax(self):
        self._ax = plt.figure().add_subplot(111)

    @property
    def _text_kwargs(self):
        return {}

    def _set_lims(self):
        self._ax.set_xlim(self.cpd.min_range, 0)
        self._ax.set_ylim(self.cpd.min_range, 0)
        self._ax.set_aspect('equal', adjustable='box')

    def _set_labels(self):
        vertex_elements = self.cpd.vertex_elements
        self._ax.set_xlabel(f"Chemical potential of {vertex_elements[0]}")
        self._ax.set_ylabel(f"Chemical potential of {vertex_elements[1]}")

    def _set_grid(self):
        self._ax.grid(**self._mpl_defaults.grid_2d)

    def draw_simplex(self, formula):
        vertex_coords = self.cpd.polygons[formula]
        plt.plot(*transpose(vertex_coords), **self._mpl_defaults.simplex_2d)


class ChemPotDiag3DMplPlotter(ChemPotDiagMplPlotter):
    def _add_ax(self):
        self._ax = plt.figure().add_subplot(111, projection='3d')

    @property
    def _text_kwargs(self):
        return {"ha": "center", "va": "center"}

    def _set_lims(self):
        self._ax.set_xlim3d(self.cpd.min_range, 0)
        self._ax.set_ylim3d(0, self.cpd.min_range)
        self._ax.set_zlim3d(self.cpd.min_range, 0)

    def _set_labels(self):
        label = [f"Chemical potential of {e}" for e in self.cpd.vertex_elements]
        self._ax.set_xlabel(label[0])
        self._ax.set_ylabel(label[1])
        self._ax.set_zlabel(label[2])

    def _set_grid(self):
        pass

    def draw_simplex(self, formula):
        vertex_coords = self.cpd.polygons[formula]
        atomic_fractions = self.cpd.atomic_fractions(formula)

        face = Poly3DCollection([sort_coords(np.array(vertex_coords))])
        face.set_color(self._3d_simplex_color(atomic_fractions))
        face.set_edgecolor("black")
        self._ax.add_collection3d(face)

    @staticmethod
    def _3d_simplex_color(atomic_fractions):
        return [(f + 1.5) / 2.5 for f in atomic_fractions]


def transpose(target_list: List[List[float]]) -> List[List[float]]:
    return np.array(target_list).transpose().tolist()


class ChemPotDiagPlotlyPlotter(ABC):
    def __init__(self, cpd: ChemPotDiag):
        self.cpd = cpd
        self.axes_titles = [f"μ <sub>{el}</sub> (eV)"
                            for el in self.cpd.vertex_elements]
        title = f"Chemical potential diagram of {self.cpd.chemical_system}"
        self.range = [self.cpd.min_range * 1.001, -self.cpd.min_range * 0.1]
        self.common_layout = dict(title=title, width=700, height=700,
                                  title_font_size=30, showlegend=False)
        self.title_axis = dict(xaxis_title=self.axes_titles[0],
                               yaxis_title=self.axes_titles[1],
                               xaxis_range=self.range,
                               yaxis_range=self.range)


class ChemPotDiag2DPlotlyPlotter(ChemPotDiagPlotlyPlotter):
    def __init__(self, cpd: ChemPotDiag):
        super().__init__(cpd)
        self.hover = "<b>label %{text}</b><br> energy (%{x:.2f}, %{y:.2f})"

    @property
    def figure(self):
        fig = go.Figure()

        for comp, vertices in self.cpd.polygons.items():
            (x0, y0), (x1, y1) = vertices[0], vertices[1]
            x_ave, y_ave = (x0 + x1) / 2, (y0 + y1) / 2
            if x_ave == 0.0:
                pos = "middle right"
            else:
                pos = "top center"

            fig.add_shape(dict(type="line", x0=x0, y0=y0, x1=x1, y1=y1,
                               line=dict(color="RoyalBlue", width=3)))
            fig.add_trace(go.Scatter(x=[x_ave], y=[y_ave],
                                     text=[clean_formula(comp)],
                                     mode="text",
                                     textposition=pos,
                                     textfont=dict(color="black", size=24),
                                     hoverinfo='skip'))

        if self.cpd.target:
            x, y, text = [], [], []
            for label, cp in self.cpd.target_coords.items():
                x.append(cp[0])
                y.append(cp[1])
                text.append(label)
            fig.add_trace(go.Scatter(x=x, y=y, text=text, mode='markers+text',
                                     textposition="bottom left",
                                     textfont=dict(size=24),
                                     hovertemplate=self.hover))

        fig.update_traces(marker_size=15)
        fig.update_layout(font_size=24, **self.common_layout, **self.title_axis)
        return fig


class ChemPotDiag3DPlotlyPlotter(ChemPotDiagPlotlyPlotter):
    def __init__(self, cpd: ChemPotDiag):
        super().__init__(cpd)
        self.hover = \
            "<b>label %{text}</b><br> energy (%{x:.2f}, %{y:.2f}, %{z:.2f})"

    @property
    def figure(self):
        data = []
        color_cycle = cycle(px.colors.qualitative.Dark24)

        for comp, vertex_coords in self.cpd.polygons.items():
            color = next(color_cycle)
            vertex_coords = sort_coords(np.array(vertex_coords))
            # add surface
            d = make_triangles(vertex_coords)
            data.append(go.Mesh3d(**d, alphahull=-1, opacity=0.3,
                                  hoverinfo='skip', color=color))
            # add formula
            data.append(go.Scatter3d(x=[sum(d["x"]) / len(vertex_coords)],
                                     y=[sum(d["y"]) / len(vertex_coords)],
                                     z=[sum(d["z"]) / len(vertex_coords)],
                                     text=[clean_formula(comp)],
                                     mode="text",
                                     textposition="middle center",
                                     textfont=dict(color=color, size=24),
                                     hoverinfo='skip'))
        if self.cpd.target:
            target_coords = self.cpd.target_coords
            x = [cp[0] for cp in target_coords.values()]
            y = [cp[1] for cp in target_coords.values()]
            z = [cp[2] for cp in target_coords.values()]
            label = [key for key in target_coords]
            data.append(go.Scatter3d(x=x, y=y, z=z, text=label,
                                     mode='markers+text',
                                     textfont=dict(size=24),
                                     hovertemplate=self.hover))
        fig = go.Figure(data=data)
        # fig.update_traces(marker_size=15)
        fig.update_layout(scene=dict(zaxis_title=self.axes_titles[2],
                                     zaxis_range=self.range,
                                     **self.title_axis),
                          font_size=15,
                          **self.common_layout)

        return fig


def clean_formula(formula: str) -> str:
    s = []
    for char in formula:
        if char == "1" or char == " ":
            continue
        elif char.isdigit():
            s.append(f"<sub>{char}</sub>")
        else:
            s.append(char)

    return "".join(s)
