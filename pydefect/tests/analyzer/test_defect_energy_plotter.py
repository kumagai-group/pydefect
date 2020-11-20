# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pydefect.analyzer.defect_energy import DefectEnergy
from pydefect.analyzer.defect_energy_plotter import DefectEnergyMplPlotter, \
    DefectEnergiesMplSettings, DefectEnergyPlotlyPlotter
from pydefect.defaults import defaults
from pymatgen.util.string import latexify
from vise.util.dash_helper import show_png


try:
    import psutil
    PSUTIL_NOT_PRESENT = False
except ModuleNotFoundError:
    PSUTIL_NOT_PRESENT = True


def test_defect_energies_mpl_settings():
    mpl_defaults = DefectEnergiesMplSettings()
    assert next(mpl_defaults.colors) == next(defaults.defect_energy_colors)

    mpl_defaults = DefectEnergiesMplSettings(colors=["black"])
    assert next(mpl_defaults.colors) == "black"


@pytest.fixture
def defect_energy_plotters():
    va_o = DefectEnergy(name="Va_O1", charges=[0, 1, 2], energies=[5, 2, -5], corrections=[1, 1, 1])
    va_mg = DefectEnergy(name="Va_Mg1", charges=[-2, -1, 0], energies=[5, 2, 0], corrections=[-1, -1, -1])
    mg_i = DefectEnergy(name="Mg_i1", charges=[1], energies=[4], corrections=[1])

    d = dict(title=latexify("MgAl2O4"),
             defect_energies=[va_o, va_mg, mg_i],
             vbm=1.5, cbm=5.5, supercell_vbm=1.0, supercell_cbm=6.0)

    return DefectEnergyMplPlotter(**d), DefectEnergyPlotlyPlotter(**d)


def test_defect_energies_mpl_actual_plot(defect_energy_plotters):
    plotter, _ = defect_energy_plotters
    plotter.construct_plot()
    plotter.plt.show()


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_defect_energies_plotly_actual_plot(defect_energy_plotters):
    _, plotter = defect_energy_plotters
    fig = plotter.create_figure()
    show_png(fig)


"""
TODO
- Show all energies with thinner lines.
- show_transition_levels

- Add labels with optimize positions
https://github.com/cphyc/matplotlib-label-lines
https://github.com/Phlya/adjustText

- Add boundary points with open circles for shallow defects.
DONE
"""