# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from matplotlib import pyplot as plt
from monty.serialization import loadfn
from pydefect.analyzer.concentration.concentration import Concentration, \
    CarrierConcentration, ConcentrationByFermiLevel
from pydefect.analyzer.concentration.plot_concentration import plot_pn, \
    DefectConcentrationMplPlotter, plot_multiple_pns


def test_plot_concentration():
    cc = [Concentration(x, CarrierConcentration(10**(-x), 10**x), [])
          for x in range(1, 10)]
    c = ConcentrationByFermiLevel(300, cc)
    plot_multiple_pns([c], [-4, 5], energy_ranges=[1, 7])

    plt.show()


def test_plot_defect_concentration(test_data_files):
    con_by_Ef = loadfn(test_data_files / "Na3AgO2" / "con_by_Ef_A_300K.json")
    from matplotlib import pyplot as plt

    plotter = DefectConcentrationMplPlotter(con_by_Ef, ax=plt.gca())
    plotter.construct_plot()
    plt.savefig("defect_concentration.pdf")
    plt.show()