# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from typing import List

from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.helpers.layouts import Column, Columns, dcc, html
from dash.dependencies import Input, Output
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import make_defect_energies
from pydefect.analyzer.defect_energy_plotter import DefectEnergyPlotlyPlotter
from pydefect.analyzer.make_defect_energy import make_single_defect_energy
from pydefect.chem_pot_diag.chem_pot_diag import CpdPlotInfo
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiagPlotly2DMplPlotter, \
    ChemPotDiagPlotly3DMplPlotter
from pydefect.corrections.abstract_correction import Correction
from pydefect.input_maker.defect_entry import DefectEntry


class CpdEnergyComponent(MPComponent):
    def __init__(self,
                 cpd_plot_info: CpdPlotInfo,
                 perfect: CalcResults,
                 defects: List[CalcResults],
                 defect_entries: List[DefectEntry],
                 corrections: List[Correction],
                 *args, **kwargs):
        self.cpd_plot_info = cpd_plot_info
        self.perfect = perfect
        self.defects = defects
        self.defect_entries = defect_entries
        self.corrections = corrections

        super().__init__(*args, **kwargs)

    @property
    def _sub_layouts(self):
        if self.cpd_plot_info.dim == 2:
            plotter = ChemPotDiagPlotly2DMplPlotter(self.cpd_plot_info)
        else:
            plotter = ChemPotDiagPlotly3DMplPlotter(self.cpd_plot_info)

        cpd = dcc.Graph(id="cpd", figure=plotter.figure)
        energy = dcc.Graph(id="energy")
        return {"cpd": cpd, "energy": energy}

    @property
    def layout(self):
        return Columns(
                    [
                        Column(html.Div(self._sub_layouts["cpd"])),
                        Column(html.Div(self._sub_layouts["energy"])),
                    ],
                    centered=True,
                )

    def generate_callbacks(self, app, cache):
        @app.callback(
            Output('energy', 'figure'),
            [Input('cpd', 'clickData')])
        def display_click_data(click_data):
            if click_data is None:
                label = "A"
            else:
                label = click_data["points"][0]["text"]
            abs_chem_pot = self.cpd_plot_info.cpd.abs_chem_pot_dict(label)

            single_energies = []
            for d, e, c in zip(self.defects, self.defect_entries, self.corrections):
                single_energies.append(
                    make_single_defect_energy(self.perfect, d, e, abs_chem_pot, c))
            defect_energies = make_defect_energies(single_energies)

            plotter = DefectEnergyPlotlyPlotter(title=f"Label {label}",
                                                defect_energies=defect_energies,
                                                vbm=self.perfect.vbm,
                                                cbm=self.perfect.cbm,
                                                supercell_vbm=self.perfect.vbm,
                                                supercell_cbm=self.perfect.cbm)

            return plotter.create_figure()

