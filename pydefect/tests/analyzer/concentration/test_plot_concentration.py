# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from matplotlib import pyplot as plt
from pydefect.analyzer.concentration.concentration import Concentration, \
    CarrierConcentration, ConcentrationByFermiLevel
from pydefect.analyzer.concentration.plot_concentration import plot_pn, \
    plot_defect_concentration


def test_plot_concentration():
    cc = [Concentration(x, CarrierConcentration(10**(-x), 10**x), [])
          for x in range(1, 10)]
    c = ConcentrationByFermiLevel(300, cc)
    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    plot_pn(c, ax1)

    plot_defect_concentration(c, ax2)

    plt.show()