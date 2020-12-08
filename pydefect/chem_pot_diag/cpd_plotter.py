# -*- coding: utf-8 -*-
from abc import abstractmethod, ABC
from itertools import cycle
from typing import Optional

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from pydefect.chem_pot_diag.chem_pot_diag import CpdPlotInfo
from pymatgen import Composition
from pymatgen.util.string import latexify
from vise.util.plotly_util import sort_coords, make_triangles


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
                 cpd_plot_info: CpdPlotInfo,
                 mpl_defaults: Optional[CpdMplSettings] = CpdMplSettings()):
        self.info = cpd_plot_info
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
        for composition, position in self.info.comp_centers.items():
            self.draw_simplex(composition)
            self._ax.text(*position,
                          latexify(composition.reduced_formula),
                          color=self._text_color(composition),
                          **self._mpl_defaults.label,
                          **self._text_kwargs)

    @abstractmethod
    def draw_simplex(self, composition):
        pass

    @property
    @abstractmethod
    def _text_kwargs(self):
        pass

    def _text_color(self, composition):
        if self.info.cpd.target and composition == self.info.cpd.target:
            return self._mpl_defaults.target_comp_text_color
        return self._mpl_defaults.comp_text_color

    def _set_title(self, title):
        if not title:
            elements_str = [str(elem) for elem in self.info.cpd.vertex_elements]
            title = f"Chemical potential diagram of {'-'.join(elements_str)}"
        self._ax.set_title(title)

    def _add_alphabet_labels(self):
        if self.info.cpd.target:
            for label, cp in self.info.cpd.target_vertices.items():
                self._ax.text(*transpose(cp), label,
                              **self._mpl_defaults.alphabet_label)


class ChemPotDiagMpl2DMplPlotter(ChemPotDiagMplPlotter):

    def _add_ax(self):
        self._ax = plt.figure().add_subplot(111)

    @property
    def _text_kwargs(self):
        return {}

    def _set_lims(self):
        self._ax.set_xlim(self.info.min_range, 0)
        self._ax.set_ylim(self.info.min_range, 0)
        self._ax.set_aspect('equal', adjustable='box')

    def _set_labels(self):
        vertex_elements = self.info.cpd.vertex_elements
        self._ax.set_xlabel(f"Chemical potential of {vertex_elements[0]}")
        self._ax.set_ylabel(f"Chemical potential of {vertex_elements[1]}")

    def _set_grid(self):
        self._ax.grid(**self._mpl_defaults.grid_2d)

    def draw_simplex(self, composition):
        vertex_coords = self.info.comp_vertices[composition]
        plt.plot(*transpose(vertex_coords), **self._mpl_defaults.simplex_2d)


class ChemPotDiagMpl3DMplPlotter(ChemPotDiagMplPlotter):
    def _add_ax(self):
        self._ax = plt.figure().add_subplot(111, projection='3d')

    @property
    def _text_kwargs(self):
        return {"ha": "center", "va": "center"}

    def _set_lims(self):
        self._ax.set_xlim3d(self.info.min_range, 0)
        self._ax.set_ylim3d(0, self.info.min_range)
        self._ax.set_zlim3d(self.info.min_range, 0)

    def _set_labels(self):
        vertex_elements = self.info.cpd.vertex_elements
        self._ax.set_xlabel(f"Chemical potential of {vertex_elements[0]}")
        self._ax.set_ylabel(f"Chemical potential of {vertex_elements[1]}")
        self._ax.set_zlabel(f"Chemical potential of {vertex_elements[2]}")

    def _set_grid(self):
        pass

    def draw_simplex(self, composition):
        vertex_coords = self.info.comp_vertices[composition]
        atomic_fractions = self.info.atomic_fractions(composition)

        face = Poly3DCollection([sort_coords(np.array(vertex_coords))])
        face.set_color(self._3d_simplex_color(atomic_fractions))
        face.set_edgecolor("black")
        self._ax.add_collection3d(face)

    @staticmethod
    def _3d_simplex_color(atomic_fractions):
        return [(f + 1.5) / 2.5 for f in atomic_fractions]


def transpose(target_list) -> list:
    return np.array(target_list).transpose().tolist()


class ChemPotDiagPlotly2DMplPlotter:

    def __init__(self, cpd_plot_info: CpdPlotInfo):
        self.info = cpd_plot_info

    @property
    def figure(self):
        fig = go.Figure()

        for comp, vertices in self.info.comp_vertices.items():
            (x0, y0), (x1, y1) = vertices[0], vertices[1]
            x_ave, y_ave = (x0 + x1) / 2, (y0 + y1) / 2
            if x_ave == 0.0:
                pos = "middle right"
            elif y_ave == 0.0:
                pos = "top center"
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

        if self.info.cpd.target:
            x, y, text = [], [], []
            for label, cp in self.info.cpd.target_vertices.items():
                x.append(cp[0])
                y.append(cp[1])
                text.append(label)
            fig.add_trace(go.Scatter(x=x, y=y, text=text, mode='markers+text',
                                     textposition="bottom left",
                                     textfont=dict(size=24),
                                     hovertemplate="<b>label %{text}</b><br>"
                                                   "energy (%{x:.2f}, %{y:.2f})"))

        vertex_elements = self.info.cpd.vertex_elements
        fig.update_traces(marker_size=15)
        elements_str = [str(elem) for elem in self.info.cpd.vertex_elements]
        fig.update_layout(
            title=f"Chemical potential diagram of {'-'.join(elements_str)}",
            xaxis_title=f"&#181;<sub>{vertex_elements[0]}</sub> (eV)",
            yaxis_title=f"&#181;<sub>{vertex_elements[1]}</sub> (eV)",
            width=700, height=700,
            title_font_size=30,
            font_size=24,
            showlegend=False)

        _range = [self.info.min_range, -self.info.min_range * 0.1]
        fig.update_xaxes(range=_range)
        fig.update_yaxes(range=_range)
        return fig


class ChemPotDiagPlotly3DMplPlotter:

    def __init__(self, cpd_plot_info: CpdPlotInfo):
        self.info = cpd_plot_info

    @property
    def figure(self):
        data = []
        color_cycle = cycle(px.colors.qualitative.Dark24)

        for comp, vertex_coords in self.info.comp_vertices.items():

            color = next(color_cycle)
            vertex_coords = sort_coords(np.array(vertex_coords))
            # add surface
            d = make_triangles(vertex_coords)
            data.append(go.Mesh3d(**d,
                                  alphahull=-1, opacity=0.3, hoverinfo='skip',
                                  color=color))

            x_ave = sum(d["x"]) / len(vertex_coords)
            y_ave = sum(d["y"]) / len(vertex_coords)
            z_ave = sum(d["z"]) / len(vertex_coords)

            # add formula
            data.append(go.Scatter3d(x=[x_ave], y=[y_ave], z=[z_ave],
                                     text=[clean_formula(comp)],
                                     mode="text",
                                     textposition="middle center",
                                     textfont=dict(color=color, size=24),
                                     hoverinfo='skip'))

        if self.info.cpd.target:
            x, y, z, text = [], [], [], []
            for label, cp in self.info.cpd.target_vertices.items():
                x.append(cp[0])
                y.append(cp[1])
                z.append(cp[2])
                text.append(label)
            data.append(go.Scatter3d(x=x, y=y, z=z, text=text,
                                     mode='markers+text',
                                     textfont=dict(size=24),
                                     hovertemplate="<b>label %{text}</b><br>"
                                                   "energy (%{x:.2f}, %{y:.2f})"))

        fig = go.Figure(data=data)
        # fig.update_traces(marker_size=15)
        _range = [self.info.min_range * 1.001, -self.info.min_range * 0.1]
        vertex_elements = self.info.cpd.vertex_elements
        elements_str = [str(elem) for elem in self.info.cpd.vertex_elements]
        fig.update_layout(
            title=f"Chemical potential diagram of {'-'.join(elements_str)}",
            scene=dict(
                xaxis_title=f"{vertex_elements[0]} (eV)",
                yaxis_title=f"{vertex_elements[1]} (eV)",
                zaxis_title=f"{vertex_elements[2]} (eV)",
                xaxis_range=_range,
                yaxis_range=_range,
                zaxis_range=_range,
            ),
            width=700, height=700,
            font_size=15,
            title_font_size=30,
            showlegend=False)

        return fig


def clean_formula(comp: Composition):
    s = []
    for char in str(comp):
        if char == "1" or char == " ":
            continue
        elif char.isdigit():
            s.append(f"<sub>{char}</sub>")
        else:
            s.append(char)

    return "".join(s)
