# coding: utf-8
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys
from pathlib import Path

import crystal_toolkit.components as ctc
from crystal_toolkit.helpers.layouts import *
from dash import Dash
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.dash_components.cpd_energy_dash import CpdEnergyComponent
from pydefect.analyzer.dash_components.scenes_from_volumetric_data import \
    SceneDicts
from pydefect.analyzer.dash_components.single_defect_component import \
    SingleDefectComponent
from pydefect.analyzer.eigenvalue_plotter import EigenvaluePlotlyPlotter
from pydefect.chem_pot_diag.chem_pot_diag import CpdPlotInfo, ChemPotDiag
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotlyPlotter
from pydefect.input_maker.defect_entry import DefectEntry

app = Dash(suppress_callback_exceptions=True,
           assets_folder=SETTINGS.ASSETS_PATH)

def create_ctk(cpd_energy_layout, single_defect_layouts):
    _layout = Container(
        [
            Section([Column(cpd_energy_layout),
                     Column(single_defect_layouts)])
        ]
    )

    return _layout


def make_layouts(perfect_dirname, defect_dirnames, chem_pot_diag):
    cpd_plot_info = CpdPlotInfo(ChemPotDiag.from_yaml(chem_pot_diag))
    perfect: CalcResults = loadfn(perfect_dirname / "calc_results.json")

    defects, defect_entries, corrections, single_defect_layouts = [], [], [], []
    for d in defect_dirnames:
        #TODO: check convergence
        defects.append(loadfn(d / "calc_results.json"))
        defect_entry = loadfn(d / "defect_entry.json")
        defect_entries.append(defect_entry)
        corrections.append(loadfn(d / "correction.json"))
        defect_entry: DefectEntry = loadfn(d / "defect_entry.json")
        efnv_correction = loadfn(d / "correction.json")
        eigval = loadfn(d / "band_edge_eigenvalues.json")

        pot = SitePotentialPlotlyPlotter(title=defect_entry.full_name,
                                         efnv_correction=efnv_correction)
        eig = EigenvaluePlotlyPlotter(title=defect_entry.full_name,
                                      band_edge_eigenvalues=eigval,
                                      supercell_vbm=perfect.vbm,
                                      supercell_cbm=perfect.cbm)
        scene_dicts = loadfn(d / "parchgs" / "scene_dicts.json")
        if isinstance(scene_dicts, dict):
            scene_dicts = SceneDicts.from_dict(scene_dicts)

        single_defect_layouts.append(
            SingleDefectComponent(pot, eig, scene_dicts, defect_entry.full_name,
                                  id=f"{defect_entry.full_name}").layout())

    cpd_energy_comp = CpdEnergyComponent(cpd_plot_info, perfect, defects,
                                         defect_entries, corrections)
    return create_ctk(cpd_energy_comp.layout, single_defect_layouts)


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-p", "--perfect", type=Path)
    parser.add_argument("-d", "--defects", type=Path, nargs="+")
    parser.add_argument("-c", "--chem_pot_diag", type=str)
    parser.add_argument("--port", type=int)
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    layout = make_layouts(args.perfect, args.defects, args.chem_pot_diag)
    ctc.register_crystal_toolkit(app=app, layout=layout, cache=None)
    app.run_server(debug=True, port=args.port)
