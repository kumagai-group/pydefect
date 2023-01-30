# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from collections import defaultdict

from matplotlib.axes import Axes
from pydefect.analyzer.concentration.concentration import \
    ConcentrationByFermiLevel


def plot_pn(cc: ConcentrationByFermiLevel, ax: Axes):
    Efs, ps, ns = [], [], []
    for cc in cc.concentrations:
        Efs.append(cc.Ef)
        ps.append(cc.carrier.p)
        ns.append(cc.carrier.n)

    ax.set_yscale("log")
    ax.plot(Efs, ps)
    ax.plot(Efs, ns)


def plot_defect_concentration(cc: ConcentrationByFermiLevel, ax: Axes):
    Efs, ddd = [], defaultdict(list)
    for cc in cc.concentrations:
        Efs.append(cc.Ef)
        for defect in cc.defects:
            ddd[defect.name].append(defect.total_concentration)

    ax.set_yscale("log")
    for dd in ddd:
        ax.plot(Efs, dd)