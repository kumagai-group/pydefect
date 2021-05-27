# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import argparse
import sys

from monty.serialization import loadfn
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag
from pydefect.chem_pot_diag.cpd_plotter import ChemPotDiag3DPlotlyPlotter, \
    ChemPotDiag2DPlotlyPlotter


def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--chem_pot_diag", dest="cpd", type=loadfn)
    parser.add_argument("--port", type=int)
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    cpd: ChemPotDiag = args.cpd
    if cpd.dim == 2:
        plotter = ChemPotDiag2DPlotlyPlotter(args.cpd)
    elif cpd.dim == 3:
        plotter = ChemPotDiag3DPlotlyPlotter(args.cpd)
    else:
        raise ValueError("Dimension must be 2 or 3.")
    fig = plotter.figure
    fig.show()
#    app.run_server(debug=True, port=args.port)