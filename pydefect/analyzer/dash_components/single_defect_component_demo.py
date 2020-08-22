# standard Dash imports
# import for this example
from pathlib import Path

# standard Crystal Toolkit import
import crystal_toolkit.components as ctc
import dash
import dash_html_components as html
from ace.dbs.utils.scenes_from_volumetric_data import SceneDicts
from ace.htmls.single_defect_component import SingleDefectComponent
from pydefect.analyzer.band_edge_states import BandEdgeEigenvalues
from pydefect.analyzer.eigenvalue_plotter import EigenvaluePlotlyPlotter
from pydefect.corrections.efnv_correction.efnv_correction import \
    ExtendedFnvCorrection, PotentialSite
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotlyPlotter
from pymatgen import loadfn

# create Dash app as normal
app = dash.Dash()

efnv_cor = ExtendedFnvCorrection(
    charge=1,
    point_charge_correction=0.0,
    defect_region_radius=2.8,
    sites=[
        PotentialSite(specie="H", distance=1.0, potential=-4, pc_potential=None),
        PotentialSite(specie="H", distance=2.0, potential=-3, pc_potential=None),
        PotentialSite(specie="He", distance=3.0, potential=-2,
                      pc_potential=-3),
        PotentialSite(specie="He", distance=4.0, potential=-1, pc_potential=-2)
    ],
    defect_coords=(0.0, 0.0, 0.0))

eig = BandEdgeEigenvalues(energies_and_occupations=[
    [[[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]], [[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]]],
    [[[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]], [[0.0, 1.0], [0.5, 0.5], [1.0, 0.0]]],
],
    kpt_coords=[(0.0, 0.0, 0.0), (0.25, 0.0, 0.0)],
    lowest_band_index=10)


potential_plotter = SitePotentialPlotlyPlotter(title="test", efnv_correction=efnv_cor)
eigenvalue_plotter = EigenvaluePlotlyPlotter(title="test", band_edge_eigenvalues=eig, supercell_vbm=3, supercell_cbm=5)

scene_dicts: SceneDicts = loadfn(str(Path(__file__).parent / ".." / "test_files" / "vasp_files" / "mgo_defects" / "actual_parchg_Va_O1_2" / "scene_dicts.json"))
#print(type(scene_dicts.structure_graph.structure))
#scene_dicts.structure_graph.structure.append(DummySpecie(), [0, 0.5, 0.5])
structure_component = SingleDefectComponent(potential_plotter, eigenvalue_plotter, scene_dicts, links={})

#structure_graph = StructureGraph.with_local_env_strategy(structure, bonding_strategy)
#structure_graph = loadfn("structure_graph.json")
#structure_component = ctc.StructureMoleculeComponent(structure_graph, id="my_structure")
#structure_component = DefectStructureComponent(structure_graph, id="my_structure")
#sg = structure_graph.to_json()
#sg_output_file = (Path(".") / "structure_graph.json")
#sg_output_file.write_text(sg)

# example layout to demonstrate capabilities of component
my_layout = html.Div(
    [
        html.H1("StructureMoleculeComponent Example"),
        html.H2("Standard Layout"),
        structure_component.layout(),
    ]
)

# tell crystal toolkit about your app and layout
ctc.register_crystal_toolkit(app, layout=my_layout)

if __name__ == "__main__":
    app.run_server(debug=True, port=8053)

