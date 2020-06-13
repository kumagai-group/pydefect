# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from monty.serialization import loadfn
from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection, DefectSite
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotter


def test_site_potential_plotter_with_simple_example():
    efnv_cor = ExtendedFnvCorrection(
        charge=1,
        point_charge_correction=0.0,
        defect_region_radius=2.8,
        sites=[
            DefectSite(specie="H", distance=1.0, potential=-4, pc_potential=-5),
            DefectSite(specie="H", distance=2.0, potential=-3, pc_potential=-4),
            DefectSite(specie="He", distance=3.0, potential=-2,
                       pc_potential=-3),
            DefectSite(specie="He", distance=4.0, potential=-1, pc_potential=-2)
                                     ])
    plotter = SitePotentialPlotter("ZnO Va_O1_2", efnv_cor)
    plotter.construct_plot()
    plotter.plt.show()


def test_site_potential_plotter_with_actual_file():
    efnv_cor = loadfn(Path(__file__).parent / "NaCl_Va_Na_-1_correction.json")
    plotter = SitePotentialPlotter("NaCl Va_Na_-1", efnv_cor)
    plotter.construct_plot()
    plotter.plt.show()

