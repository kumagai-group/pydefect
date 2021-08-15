# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import dash_table
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.helpers.layouts import Column, dcc, html, Columns
from dash.dependencies import Input, Output
from pydefect.analyzer.defect_energy import DefectEnergySummary
from pydefect.analyzer.defect_energy_plotter import DefectEnergyPlotlyPlotter
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, TargetVertices
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiag2DPlotlyPlotter, \
    ChemPotDiag3DPlotlyPlotter


def _target_vertices_to_columns_and_data(chem_pot_diag: ChemPotDiag):
    chem_pots = chem_pot_diag.to_target_vertices.chem_pots
    columns = [{"name": "label", "id": "label"}]
    for elem in chem_pot_diag.vertex_elements:
        columns.append({"name": elem, "id": elem})

    data = []
    for k, v in chem_pots.items():
        datum = {"label": k}
        for elem in chem_pot_diag.vertex_elements:
            datum[elem] = v[elem]
        data.append(datum)
    return columns, data


class CpdEnergyComponent(MPComponent):
    def __init__(self,
                 chem_pot_diag: ChemPotDiag,
                 defect_energy_summary: DefectEnergySummary,
                 *args, **kwargs):
        self.chem_pot_diag = chem_pot_diag
        self.dim = chem_pot_diag.dim
        self.defect_energy_summary = defect_energy_summary
        super().__init__(*args, **kwargs)

        columns, data = _target_vertices_to_columns_and_data(chem_pot_diag)
        self.vertex_list = dash_table.DataTable(
            id=self.id("vertex_list"),
            columns=columns,
            data=data,
            style_cell={'fontSize': 20, 'font-family': 'sans-serif',
                        "width": 700})

    def energy_fig(self, label, allow_shallow):
        plotter = DefectEnergyPlotlyPlotter(
            title=f"Label {label[0]}",
            chem_pot_label=label,
            defect_energy_summary=self.defect_energy_summary,
            allow_shallow=allow_shallow,
            with_corrections=True)
        return plotter.create_figure()

    @property
    def _sub_layouts(self):
        allow_shallow = self.get_choice_input(
            kwarg_label="allow_shallow",
            state={"label": "True", "value": "True"},
            label="Allow defects with perturbed host states",
            # Only string values are allowed.
            options=[{"label": "True", "value": "True"},
                     {"label": "False", "value": "False"}])

        if self.dim in [2, 3]:
            if self.dim == 2:
                fig = ChemPotDiag2DPlotlyPlotter(self.chem_pot_diag).figure
            else:
                fig = ChemPotDiag3DPlotlyPlotter(self.chem_pot_diag).figure
            cpd = dcc.Graph(id="cpd", figure=fig)
            energy = dcc.Graph(id=self.id('energy'))
            return {"cpd": cpd, "energy": energy, "allow_shallow": allow_shallow}
        else:
            d = [{"label": v, "value": v}
                 for v in self.chem_pot_diag.target_vertices_dict.keys()]
            cpd_label = self.get_choice_input(
                kwarg_label="label",
                state={"label": "A", "value": "A"},
                label="Equilibrium label",
                options=d)
            energy = dcc.Graph(id=self.id('energy'))
            return {"cpd_label": cpd_label, "energy": energy,
                    "allow_shallow": allow_shallow}

    @property
    def layout(self):
        if self.dim in [2, 3]:
            return Columns([Column([html.Div(self._sub_layouts["energy"]),
                                    self._sub_layouts["allow_shallow"]]),
                            Column([self._sub_layouts["cpd"],
                                    self.vertex_list])])

        else:
            return Columns([Column([self.vertex_list,
                                    self._sub_layouts["allow_shallow"],
                                    self._sub_layouts["cpd_label"]]),
                            Column(html.Div(self._sub_layouts["energy"]))])

    def generate_callbacks(self, app, cache):
        if self.dim in [2, 3]:
            @app.callback(
                Output(self.id('energy'), 'figure'),
                [Input('cpd', 'clickData'),
                 Input(self.get_kwarg_id('allow_shallow'), 'value')])
            def display_cpd(click_data, allow_shallow):
                allow_shallow = True if allow_shallow[0] == "True" else False
                label = click_data["points"][0]["text"] if click_data else "A"
                return self.energy_fig(label, allow_shallow)
        else:
            @app.callback(
                Output(self.id('energy'), 'figure'),
                [Input(self.get_kwarg_id("label"), "value"),
                 Input(self.get_kwarg_id('allow_shallow'), 'value')])
            def display_cpd(label, allow_shallow):
                allow_shallow = True if allow_shallow[0] == "True" else False
                label = label[0] if label else "A"
                return self.energy_fig(label, allow_shallow)


