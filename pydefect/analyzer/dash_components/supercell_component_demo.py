# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

# standard Dash imports
# import for this example
import argparse
import sys

# standard Crystal Toolkit import
import crystal_toolkit.components as ctc
import dash
import dash_html_components as html
from pydefect.analyzer.dash_components.supercell_component import \
    SupercellComponent
from pymatgen import loadfn


# example layout to demonstrate capabilities of component
def parse_args(args):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-s", "--supercell_info", type=loadfn)
    return parser.parse_args(args)


if __name__ == "__main__":
    app = dash.Dash()
    args = parse_args(sys.argv[1:])
    supercell_component = SupercellComponent(args.supercell_info)
    my_layout = html.Div(
        [
            html.H1("SupercellComponent Example"),
            html.H2("Standard Layout"),
            supercell_component.layout,
        ]
    )

    ctc.register_crystal_toolkit(app, layout=my_layout)
    app.run_server(debug=True, port=8053)



