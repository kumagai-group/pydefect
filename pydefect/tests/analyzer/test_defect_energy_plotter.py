# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest
from pydefect.analyzer.defect_energy import DefectEnergy, DefectEnergies, \
    DefectEnergySummary
from pydefect.analyzer.defect_energy_plotter import DefectEnergyMplPlotter, \
    DefectEnergiesMplSettings, DefectEnergyPlotlyPlotter
from pydefect.defaults import defaults
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
def defect_energy_summary():
    de1 = DefectEnergy(0.0, {"PC correction": 2.0}, False)
    de2 = DefectEnergy(-5.0, {"PC correction": 3.0}, True)

    defect_energies = {"Va_Mg1": DefectEnergies(atom_io={"Mg": -1},
                                                charges=[0, 1],
                                                defect_energies=[de1, de2])}
    return DefectEnergySummary(title="MgAl2O4",
                               defect_energies=defect_energies,
                               rel_chem_pots={"A": {"Mg": -1}},
                               cbm=7.0,
                               supercell_vbm=-1.0,
                               supercell_cbm=4.0)


def test_defect_energies_mpl_plot(defect_energy_summary):
    plotter = DefectEnergyMplPlotter(
        defect_energy_summary=defect_energy_summary,
        chem_pot_label="A",
        allow_shallow=True,
        with_corrections=True,
        add_thin_lines=True)
    plotter.construct_plot()
    plotter.plt.show()


@pytest.mark.skipif(PSUTIL_NOT_PRESENT, reason="psutil does not exist")
def test_defect_energies_plotly_plot(defect_energy_summary):
    plotter = DefectEnergyPlotlyPlotter(
        defect_energy_summary=defect_energy_summary,
        chem_pot_label="A",
        allow_shallow=False,
        with_corrections=True)
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