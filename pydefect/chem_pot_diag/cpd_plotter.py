# -*- coding: utf-8 -*-
from abc import abstractmethod, ABC
from typing import Optional

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numpy import concatenate, clip, dot, arctan2
from numpy.linalg import det
from pydefect.chem_pot_diag.chem_pot_diag import CpdPlotInfo
from pymatgen.util.string import latexify


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


class ChemPotDiagPlotter(ABC):
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


def transpose(target_list) -> list:
    return np.array(target_list).transpose().tolist()


class ChemPotDiag2DPlotter(ChemPotDiagPlotter):

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


class ChemPotDiag3DPlotter(ChemPotDiagPlotter):
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


def sort_coords(coords: np.ndarray) -> np.ndarray:
    """Sort coordinates based on the angle with first coord from the center.

    Args:
        coords (np.ndarray):
            Coordinates to be sorted. The format of coords is as follows.
            np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]

    Returns:
        np.ndarray for sorted coordinates.
    """
    if len(coords[0]) != 3:
        raise ValueError("Only valid for 3D vector")

    center = np.average(coords, axis=0)
    relative_coords = coords - center
    external_prod = np.cross(relative_coords[0], relative_coords[1])
    if abs(np.linalg.norm(external_prod)) < 1e-8:  # Skip parallel vectors.
        external_prod = np.cross(relative_coords[0], relative_coords[2])
    normal_to_12_plane = external_prod / np.linalg.norm(external_prod)

    v0 = relative_coords[0] / np.linalg.norm(relative_coords[0])

    def angle_between_v0(index: int) -> float:
        """
        Args:
            index (int): index of coords.

        Returns (float):
            Angle between rays from the center to rel_coords[0] and
            rel_coords[int].
        """
        v = relative_coords[index] / np.linalg.norm(relative_coords[index])
        matrix = concatenate(([v0], [v], [normal_to_12_plane]), axis=0)
        determinant = det(matrix)
        angle = arctan2(clip(dot(v0, v), -1.0, 1.0), determinant)
        return angle

    indices = [i for i in range(len(coords))]
    indices.sort(key=angle_between_v0)
    return coords[indices]
