# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import pytest
from monty.serialization import loadfn
from pydefect.corrections.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.corrections.site_potential_plotter import \
    SitePotentialMplPlotter, SitePotentialPlotlyPlotter
from vise.util.dash_helper import show_png

try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


@pytest.fixture
def efnv_cor():
    return ExtendedFnvCorrection(
        charge=1,
        point_charge_correction=0.0,
        defect_region_radius=2.8,
        sites=[
            PotentialSite(specie="H", distance=1.0, potential=-4, pc_potential=None),
            PotentialSite(specie="H", distance=2.0, potential=-3, pc_potential=None),
            PotentialSite(specie="He", distance=3.0, potential=-2,
                          pc_potential=-3),
            PotentialSite(specie="He", distance=4.0, potential=-1, pc_potential=-2)
        ],
        defect_coords=(0.0, 0.0, 0.0))


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="skipped for circle CI")
def test_site_potential_plotter_with_simple_example(efnv_cor):
    plotter = SitePotentialMplPlotter.from_efnv_corr(
        title="ZnO Va_O1_2", efnv_correction=efnv_cor)
    plotter.construct_plot()
    plotter.plt.show()


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="skipped for circle CI")
def test_site_potential_plotter_with_actual_file():
    efnv_cor = loadfn(Path(__file__).parent / "NaCl_Va_Na_-1_correction.json")
    plotter = SitePotentialMplPlotter.from_efnv_corr(
        title="NaCl Va_Na_-1", efnv_correction=efnv_cor)
    plotter.construct_plot()
    plotter.plt.show()


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_defect_energies_plotly_actual_plot(efnv_cor):
    fig = SitePotentialPlotlyPlotter.from_efnv_corr(
        title="ZnO Va_O1_2", efnv_correction=efnv_cor).create_figure()
    show_png(fig)


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="skipped for circle CI")
def test_gkfo_site_potential_plotter():
    gkfo = loadfn(Path(__file__).parent / "gkfo_correction.json")
    plotter = SitePotentialMplPlotter.from_gkfo_corr(
        title="GKFO", gkfo_correction=gkfo)
    plotter.construct_plot()
    plotter.plt.show()


