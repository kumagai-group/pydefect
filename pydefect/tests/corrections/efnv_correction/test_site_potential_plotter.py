# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import pytest

from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection, DefectSite
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotter


@pytest.fixture
def plotter():
    efnv_cor = ExtendedFnvCorrection(charge=1,
                                 point_charge_correction=0.0,
                                 defect_region_radius=2.8,
                                 sites=[
                                     DefectSite(specie="H", distance=1.0, potential=-4, pc_potential=-5),
                                     DefectSite(specie="H", distance=2.0, potential=-3, pc_potential=-4),
                                     DefectSite(specie="He", distance=3.0, potential=-2, pc_potential=-3),
                                     DefectSite(specie="He", distance=4.0, potential=-1, pc_potential=-2)
                                 ]
                                 )
    return SitePotentialPlotter("ZnO Va_O1_2", efnv_cor)


def test(plotter):
    plotter.construct_plot()
    plotter.plt.show()


"""
TODO
- Plot dft potentials

DONE
"""