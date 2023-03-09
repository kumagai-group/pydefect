# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from collections import defaultdict
from typing import List

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from pydefect.analyzer.concentration.concentration import \
    ConcentrationByFermiLevel
from vise.util.matplotlib import float_to_int_formatter


def plot_multiple_pns(ccs: List[ConcentrationByFermiLevel],
                      concentration_ranges: List[float] = None,
                      energy_ranges: List[float] = None):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    for cc, style in zip(ccs, ["-", "-.", "--", ":"]):
        plot_pn(cc, ax, style)

    if concentration_ranges:
        ax.set_ylim([10 ** i for i in concentration_ranges])

    if energy_ranges:
        ax.set_xlim(energy_ranges)

    plt.xlabel("Fermi level (eV)")
    plt.ylabel("Concentration")
    plt.legend()


def plot_pn(cc: ConcentrationByFermiLevel, ax: Axes, style: str = "-"):
    Efs, ps, ns = [], [], []
    for c in cc.concentrations:
        Efs.append(c.Ef)
        ps.append(c.carrier.p)
        ns.append(c.carrier.n)

    ax.set_yscale("log")
    ax.plot(Efs, ps, color="red", linestyle=style, label=cc.T)
    ax.plot(Efs, ns, color="blue", linestyle=style)
    ax.legend()


class DefectConcentrationMplPlotter:
    def __init__(self,
                 con_by_Ef: ConcentrationByFermiLevel,
                 **plot_settings):
        self.con_by_Ef = con_by_Ef
        self.plt = plt
        self.ax = plt.gca()

    def construct_plot(self):
        self._plot_defect_concentration()
        self._set_scale()
        self._set_labels()
        self._set_legend()
        self._set_formatter()

    def _plot_defect_concentration(self):
        Efs, ddd = [], defaultdict(list)
        for cc in self.con_by_Ef.concentrations:
            Efs.append(cc.Ef)
            for defect in cc.defects:
                ddd[defect.name].append(defect.total_concentration)

        for name, dd in ddd.items():
            self.ax.plot(Efs, dd, label=name)

    def _set_scale(self):
        self.ax.set_yscale("log")

    def _set_labels(self):
        self.plt.xlabel(f"Fermi level (eV)")
        self.plt.ylabel("Concentration (cm$^{-3}$)")

    def _set_formatter(self):
        axis = self.plt.gca()
        axis.xaxis.set_major_formatter(float_to_int_formatter)

    def _set_legend(self):
        self.ax.legend(bbox_to_anchor=(0.9, 0.5), loc='center left')