# standard Dash imports
# import for this example
import argparse
import sys
from pathlib import Path

# standard Crystal Toolkit import
import crystal_toolkit.components as ctc
import dash
import dash_html_components as html
from ace.htmls.single_defect_component import SingleDefectComponent
from monty.serialization import loadfn
from pydefect.analyzer.dash_components.scenes_from_volumetric_data import \
    SceneDicts
from pydefect.analyzer.eigenvalue_plotter import EigenvaluePlotlyPlotter
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotlyPlotter


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-p", "--perfect", type=Path, required=True)
    parser.add_argument("-d", "--defect", type=Path, required=True)
    return parser.parse_args(args)


if __name__ == "__main__":
    app = dash.Dash()
    args = parse_args(sys.argv[1:])
    efnv_cor = loadfn(str(args.defect / "correction.json"))
    eig = loadfn(args.defect / "band_edge_eigenvalues.json")
    scene_dicts = loadfn(args.defect / "parchgs" / "scene_dicts.json")
    if isinstance(scene_dicts, dict):
        scene_dicts = SceneDicts.from_dict(scene_dicts)

    defect_calc_results = loadfn(args.perfect / "calc_results.json")
    perfect = loadfn(args.perfect / "calc_results.json")

    potential_plotter = SitePotentialPlotlyPlotter(title=str(args.defect),
                                                   efnv_correction=efnv_cor)
    eigenvalue_plotter = EigenvaluePlotlyPlotter(title=str(args.defect),
                                                 band_edge_eigenvalues=eig,
                                                 supercell_vbm=perfect.vbm,
                                                 supercell_cbm=perfect.cbm)
    structure_component = SingleDefectComponent(potential_plotter, eigenvalue_plotter, scene_dicts, links={})

    my_layout = html.Div(
        [
            html.H1("StructureMoleculeComponent Example"),
            html.H2("Standard Layout"),
            structure_component.layout(),
        ]
    )

    # tell crystal toolkit about your app and layout
    ctc.register_crystal_toolkit(app, layout=my_layout)
    app.run_server(debug=True, port=8055)


