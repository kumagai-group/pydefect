#!/usr/bin/env python
# coding: utf-8

# In[1]:


from crystal_toolkit.helpers.layouts import Columns, Column
from crystal_toolkit.settings import SETTINGS
from jupyter_dash import JupyterDash
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.dash_components.cpd_energy_dash import CpdEnergyComponent
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo,     CompositionEnergy
from pydefect.corrections.manual_correction import ManualCorrection
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Composition, Structure, Lattice, Element
import dash_html_components as html
import crystal_toolkit.components as ctc
from dash.dependencies import Input, Output, State
import json

app = JupyterDash(suppress_callback_exceptions=True,
               assets_folder=SETTINGS.ASSETS_PATH)
from vise.analyzer.band_edge_properties import BandEdge

comp_energies = [
    CompositionEnergy(Composition("Mg"), 0.0, "a"),
    CompositionEnergy(Composition("Ca"), 0.0, "a"),
    CompositionEnergy(Composition("Sr"), 0.0, "a"),
    CompositionEnergy(Composition("O"), 0.0, "a"),
    CompositionEnergy(Composition("H"), 0.0, "a"),
#    CompositionEnergy(Composition("MgCaO3"), -100.0, "a"),
    CompositionEnergy(Composition("MgCaSrO3"), -100.0, "a"),
]
#cpd = ChemPotDiag(comp_energies, target=Composition("MgCaO3"))
cpd = ChemPotDiag(comp_energies, target=Composition("MgCaSrO3"))
cpd_plot_info = CpdPlotInfo(cpd)


# In[2]:


print(cpd.target.elements)
print(cpd.dim)
print(cpd.target_vertices)
print(cpd.all_compounds)
print(cpd.impurity_abs_energy(Element.H, label="A"))


# In[3]:


structure = Structure(Lattice.cubic(1), species=["O"] * 2, coords=[[0]*3]*2)
defect_structure = Structure(Lattice.cubic(1), species=["O"] * 1, coords=[[0]*3])

common = dict(site_symmetry="1",
             magnetization=0.0,
             kpoint_coords=[[0]*3],
             kpoint_weights=[1.0],
             potentials=[0.0],
             vbm_info=BandEdge(0.0),
             cbm_info=BandEdge(1.0),
             fermi_level=0.0)

perfect = CalcResults(structure=structure,energy=0, **common)
defects = [CalcResults(structure=defect_structure, energy=1.0, **common),
          CalcResults(structure=defect_structure, energy=0.5, **common)]

de_common = dict(name="Va_O1",
                structure=defect_structure, site_symmetry="1",
                perturbed_structure=defect_structure, defect_center=[[0]*3])

defect_entries = [DefectEntry(charge=0, **de_common),
                 DefectEntry(charge=1, **de_common)]

corrections = [ManualCorrection(correction_energy=1.0),
              ManualCorrection(correction_energy=1.0)]

cpd_e_component = CpdEnergyComponent(cpd_plot_info,
                                     perfect,
                                     defects,
                                     defect_entries,
                                     corrections,
                                     unitcell_vbm=0.0,
                                     unitcell_cbm=1.0)

my_layout = html.Div([Column(cpd_e_component.layout)])
ctc.register_crystal_toolkit(app=app, layout=my_layout, cache=None)


# In[5]:


app.run_server(port=8099)
#app.run_server(mode='inline', port=8096)

