#!/usr/bin/env python
# coding: utf-8

import crystal_toolkit.components as ctc
import dash_html_components as html
from crystal_toolkit.settings import SETTINGS
from jupyter_dash import JupyterDash
from pydefect.analyzer.dash_components.cpd_energy_dash import CpdEnergyComponent
from pydefect.analyzer.defect_energy import DefectEnergySummary, DefectEnergies, \
    DefectEnergy
from pydefect.chem_pot_diag.chem_pot_diag import RelativeEnergies, \
    ChemPotDiagMaker


if __name__ == "__main__":
    app = JupyterDash(suppress_callback_exceptions=True,
                      assets_folder=SETTINGS.ASSETS_PATH)

#    relative_energies = RelativeEnergies({"MgCaO2": -100.0})
    relative_energies = RelativeEnergies({"MgCaBaO2": -100.0})
    # cpd = ChemPotDiagMaker(relative_energies,
    #                        elements=["Mg", "Ca", "O"],
    #                        target="MgCaO2").chem_pot_diag
    cpd = ChemPotDiagMaker(relative_energies,
                           elements=["Mg", "Ca", "Ba", "O"],
                           target="MgCaBaO2").chem_pot_diag

    energy = DefectEnergy(0.0, {"corr": 1.0}, is_shallow=False)
    defect_energy_summary = DefectEnergySummary(
        title="test",
        defect_energies={"Va_O1": DefectEnergies(atom_io={"O": -1},
                                                 charges=[0],
                                                 defect_energies=[energy])},
        rel_chem_pots=cpd.to_target_vertices.chem_pots,
        cbm=2.0, supercell_vbm=-0.2, supercell_cbm=2.2 )

    cpd_e_component = CpdEnergyComponent(cpd, defect_energy_summary)

    my_layout = html.Div([cpd_e_component.layout])
    ctc.register_crystal_toolkit(app=app, layout=my_layout, cache=None)

    app.run_server(port=8099)

