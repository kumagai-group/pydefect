# coding: utf-8
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys
from pathlib import Path

import crystal_toolkit.components as ctc
from ace.htmls.components.supercell_component import \
    SupercellComponent
from crystal_toolkit.helpers.layouts import *
from dash import Dash
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.dash_components.cpd_energy_dash import \
    CpdEnergy2D3DComponent, CpdEnergyOtherComponent
from pydefect.analyzer.dash_components.scenes_from_volumetric_data import \
    SceneDicts
from pydefect.analyzer.dash_components.single_defect_component import \
    SingleDefectComponent
from pydefect.analyzer.eigenvalue_plotter import EigenvaluePlotlyPlotter
from pydefect.chem_pot_diag.chem_pot_diag import CpdPlotInfo, ChemPotDiag
from pydefect.corrections.site_potential_plotter import \
    SitePotentialPlotlyPlotter
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen.core import Structure
from vise.analyzer.dash_components.band_dos_dash import BandDosComponent
from vise.analyzer.dash_components.main import symmetry_layout, site_layout, \
    mpid_and_link
from vise.analyzer.dash_components.structure_component import StructureComponent
from vise.util.structure_symmetrizer import StructureSymmetrizer

app = Dash(suppress_callback_exceptions=True,
           assets_folder=SETTINGS.ASSETS_PATH)


def create_ctk(struct_component,
               crystal_symmetry,
               mpid_and_link,
               sites,
               band_dos_component,
               supercell_component_layout,
               cpd_energy_layout,
               single_defect_layouts):
    box_size = "30vmin"
    supercell = Column(
        [
            Box(
                supercell_component_layout,
                style={
                    "width": box_size,
                    "height": box_size,
                    "minWidth": "300px",
                    "minHeight": "300px",
                    "maxWidth": "600px",
                    "maxHeight": "600px",
                    "overflow": "hidden",
                    "padding": "0.25rem",
                    "marginBottom": "0.5rem",
                },
            ),
        ],
        narrow=True,
    )


    _layout = Container(
        [
            Section(
                [
                    Columns(
                        [
                            Column(
                                [struct_component.title_layout()]
                            )
                        ]
                    ),
                    Columns(
                        [
                            Column(
                                [
                                    Box(
                                        struct_component.layout(size="100%"),
                                        style={
                                            "width": box_size,
                                            "height": box_size,
                                            "minWidth": "300px",
                                            "minHeight": "300px",
                                            "maxWidth": "800px",
                                            "maxHeight": "800px",
                                            "overflow": "hidden",
                                            "padding": "0.25rem",
                                            "marginBottom": "0.5rem",
                                        },
                                    ),
                                    html.Div(
                                        [
                                            html.Div(
                                                struct_component.legend_layout(),
                                                style={"float": "left"},
                                            ),
                                        ],
                                        style={
                                            "width": box_size,
                                            "minWidth": "300px",
                                            "marginBottom": "40px",
                                        },
                                    ),
                                ],
                                narrow=True,
                            ),
                            Column(
                                [crystal_symmetry, mpid_and_link],
                                style={"width": box_size, "max-width": box_size},
                            ),
                        ],
                        desktop_only=False,
                        centered=False,
                    ),
                    Columns(Column([sites])),
                ]
            ),
            Section([band_dos_component.layout()]
                    ),
            Section([Columns(cpd_energy_layout + [supercell]),
                     Column(single_defect_layouts)])
        ]
    )

    return _layout


def make_layouts(structure: Structure, dos_plot_data, band_plot_data,
                 perfect_dirname, defect_dirnames, supercell_info, chem_pot_diag):
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

    if cpd_plot_info.cpd.dim in [2, 3]:
        cpd_energy_comp = CpdEnergy2D3DComponent(cpd_plot_info, perfect, defects,
                                             defect_entries, corrections)
    else:
        cpd_energy_comp = CpdEnergyOtherComponent(cpd_plot_info, perfect, defects,
                                                  defect_entries, corrections)

    structure_component = StructureComponent(structure)
    comp = structure.composition.reduced_formula
    band_dos_component = BandDosComponent(dos_plot_data, band_plot_data, id=f"band_dos_{comp}",)
    supercell_component = SupercellComponent(supercell_info)
    symmetrizer = StructureSymmetrizer(structure)
    return create_ctk(structure_component,
                      symmetry_layout(structure),
                      mpid_and_link(symmetrizer),
                      site_layout(symmetrizer),
                      band_dos_component,
                      supercell_component.layout,
                      cpd_energy_comp.layout,
                      single_defect_layouts)


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-b", "--band", type=Path, required=True)
    parser.add_argument("-d", "--dos", type=Path, required=True)
    parser.add_argument("-p", "--perfect", type=Path, required=True)
    parser.add_argument("-def", "--defects", type=Path, nargs="+", required=True)
    parser.add_argument("-s", "--supercell_info", type=loadfn, required=True)
    parser.add_argument("-c", "--chem_pot_diag", type=str, required=True)
    parser.add_argument("--port", type=int)
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    structure = Structure.from_file(args.dos / defaults.contcar)
    dos_plot_data = loadfn(args.dos / "dos_plot_data.json")
    band_plot_data = loadfn(args.band / "band_plot_info.json")
    layout = make_layouts(structure, dos_plot_data, band_plot_data,
                          args.perfect, args.defects, args.supercell_info,
                          args.chem_pot_diag)
    ctc.register_crystal_toolkit(app=app, layout=layout, cache=None)
    app.run_server(debug=True, port=args.port)
