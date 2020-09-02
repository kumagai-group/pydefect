# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys

from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiagPlotly3DMplPlotter


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--chem_pot_diag", type=str)
    parser.add_argument("--port", type=int)
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    cpd = ChemPotDiag.from_yaml(args.chem_pot_diag)
    print(cpd.target_vertices)
    print(cpd.vertex_coords)
    cpd_plot_info = CpdPlotInfo(cpd)
    print(cpd_plot_info.comp_vertices)
    plotter = ChemPotDiagPlotly3DMplPlotter(cpd_plot_info)
    fig = plotter.figure
    fig.show()
#    ctc.register_crystal_toolkit(app=app, layout=layout, cache=None)
#    app.run_server(debug=True, port=args.port)