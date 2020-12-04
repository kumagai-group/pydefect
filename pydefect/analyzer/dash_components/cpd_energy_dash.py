# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from typing import List

import dash_table
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.helpers.layouts import Column, dcc, html, Columns
from dash.dependencies import Input, Output
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import make_energies
from pydefect.analyzer.defect_energy_plotter import DefectEnergyPlotlyPlotter
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
                 unitcell_vbm: float,
                 unitcell_cbm: float,
                 band_edge_states: List[BandEdgeStates] = None,
                 *args, **kwargs):
        self.cpd_plot_info = cpd_plot_info
        self.perfect = perfect
        self.defects = defects
        self.defect_entries = defect_entries
        self.corrections = corrections
        self.unitcell_vbm = unitcell_vbm
        self.unitcell_cbm = unitcell_cbm
        self.band_edge_states = band_edge_states
        self.dim = cpd_plot_info.dim
        super().__init__(*args, **kwargs)

        df = self.cpd_plot_info.cpd.target_vertex_list_dataframe.copy()
        df.reset_index(inplace=True)

        def rename(x):
            if "mu" in x:
                return x.split("_")[1]
            return x

        self.vertex_list = dash_table.DataTable(
            id='table',
            columns=[{"name": rename(i), "id": i} for i in df.columns],
            data=df.to_dict('records'),
            style_cell={'fontSize': 24,
                        'font-family': 'sans-serif',
                        "width": 700})

        self.allow_shallow = self.get_choice_input(
            kwarg_label="allow_shallow",
            state={"label": "True", "value": True},
            label="Allow shallow states",
            options=[{"label": "True", "value": True},
                     {"label": "False", "value": False}],
            style={"fontSize": 18, 'font-family': 'sans-serif'})

    def energy_fig(self, abs_chem_pot, label, allow_shallow):
        defect_energies = make_energies(self.perfect,
                                        self.defects,
                                        self.defect_entries,
                                        self.corrections,
                                        abs_chem_pot,
                                        allow_shallow,
                                        self.band_edge_states)
        plotter = DefectEnergyPlotlyPlotter(title=f"Label {label[0]}",
                                            defect_energies=defect_energies,
                                            vbm=self.unitcell_vbm,
                                            cbm=self.unitcell_cbm,
                                            supercell_vbm=self.perfect.vbm,
                                            supercell_cbm=self.perfect.cbm,
                                            supercell_edge=True)
        return plotter.create_figure()

    @property
    def _sub_layouts(self):
        if self.dim in [2, 3]:
            if self.cpd_plot_info.dim == 2:
                fig = ChemPotDiagPlotly2DMplPlotter(self.cpd_plot_info).figure
            else:
                fig = ChemPotDiagPlotly3DMplPlotter(self.cpd_plot_info).figure
            cpd = dcc.Graph(id="cpd", figure=fig)
            energy = dcc.Graph(id="energy")
            return {"cpd": cpd, "energy": energy}
        else:
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
        if self.dim in [2, 3]:
            return Columns([Column([html.Div(self._sub_layouts["energy"]),
                                    self.allow_shallow]),
                            Column([self._sub_layouts["cpd"],
                                    self.vertex_list])
                            ])

        else:
            return Columns([Column([self.vertex_list,
                                    self._sub_layouts["cpd_label"]]),
                            Column(html.Div(self._sub_layouts["energy"]))])

    def generate_callbacks(self, app, cache):
        if self.dim in [2, 3]:
            @app.callback(
                Output('energy', 'figure'),
                [Input('cpd', 'clickData'),
                 Input(self.get_kwarg_id('allow_shallow'), 'value')])
            def display_cpd(click_data, allow_shallow):
                label = click_data["points"][0]["text"] if click_data else "A"
                abs_chem_pot = self.cpd_plot_info.cpd.abs_chem_pot_dict(label)
                return self.energy_fig(abs_chem_pot, label, allow_shallow[0])
        else:
            @app.callback(
                Output('energy', 'figure'),
                [Input(self.get_kwarg_id("label"), "value"),
                 Input(self.get_kwarg_id('allow_shallow'), 'value')])
            def display_cpd(label, allow_shallow):
                label = label if label else "A"
                abs_chem_pot = self.cpd_plot_info.cpd.abs_chem_pot_dict(label)
                return self.energy_fig(abs_chem_pot, label, allow_shallow[0])


