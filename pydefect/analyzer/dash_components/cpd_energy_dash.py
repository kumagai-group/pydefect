# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from typing import List

import dash_table
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.helpers.layouts import Column, dcc, html, Columns
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
        df = self.cpd_plot_info.cpd.target_vertex_list_dataframe.copy()
        df.reset_index(inplace=True)
        self.vertex_list = dash_table.DataTable(
            id='table',
            columns=[{"name": i, "id": i} for i in df.columns],
            data=df.to_dict('records'))

        super().__init__(*args, **kwargs)

    def energy_fig(self, abs_chem_pot, label):
        single_energies = []
        for d, e, c in zip(self.defects,
                           self.defect_entries,
                           self.corrections):
            single_energies.append(
                make_single_defect_energy(self.perfect, d, e, abs_chem_pot,
                                          c))
        defect_energies = make_defect_energies(single_energies)
        plotter = DefectEnergyPlotlyPlotter(title=f"Label {label[0]}",
                                            defect_energies=defect_energies,
                                            vbm=self.perfect.vbm,
                                            cbm=self.perfect.cbm,
                                            supercell_vbm=self.perfect.vbm,
                                            supercell_cbm=self.perfect.cbm)
        return plotter.create_figure()


class CpdEnergyOtherComponent(CpdEnergyComponent):
    @property
    def _sub_layouts(self):
        d = [{"label": v, "value": v}
             for v in self.cpd_plot_info.cpd.target_vertices.keys()]
        cpd_label = self.get_choice_input(
            kwarg_label="label",
            state={"label": "A", "value": "A"},
            label="Equilibrium label",
            options=d)
        energy = dcc.Graph(id="energy")
        return {"cpd_label": cpd_label, "energy": energy}

    @property
    def layout(self):
        return Columns([Column([self.vertex_list,
                                self._sub_layouts["cpd_label"]]),
                        Column(html.Div(self._sub_layouts["energy"]))])

    def generate_callbacks(self, app, cache):
        @app.callback(
            Output('energy', 'figure'),
            [Input(self.get_kwarg_id("label"), "value")])
        def display_cpd(label):
            label = label if label else "A"
            abs_chem_pot = self.cpd_plot_info.cpd.abs_chem_pot_dict(label[0])
            return self.energy_fig(abs_chem_pot, label)


class CpdEnergy2D3DComponent(CpdEnergyComponent):
    @property
    def _sub_layouts(self):

        if self.cpd_plot_info.dim == 2:
            fig = ChemPotDiagPlotly2DMplPlotter(self.cpd_plot_info).figure
        else:
            fig = ChemPotDiagPlotly3DMplPlotter(self.cpd_plot_info).figure
        cpd = dcc.Graph(id="cpd", figure=fig)
        energy = dcc.Graph(id="energy")
        return {"cpd": cpd, "energy": energy}

    @property
    def layout(self):
        return [Column([self.vertex_list, self._sub_layouts["cpd"]]),
                Column(html.Div(self._sub_layouts["energy"]))]

    def generate_callbacks(self, app, cache):
        @app.callback(
            Output('energy', 'figure'),
            [Input('cpd', 'clickData')])
        def display_cpd(click_data):
            label = click_data["points"][0]["text"] if click_data else "A"
            abs_chem_pot = self.cpd_plot_info.cpd.abs_chem_pot_dict(label)
            return self.energy_fig(abs_chem_pot, label)

