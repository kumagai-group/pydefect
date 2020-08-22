# coding: utf-8
#  Copyright (c) 2020 Kumagai group.

from typing import Optional, Tuple

from ace.dbs.utils.scenes_from_volumetric_data import SceneDicts
from crystal_toolkit.core.legend import Legend
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.core.scene import Scene
from crystal_toolkit.helpers.layouts import *
from crystal_toolkit.renderables.structuregraph import StructureGraph, \
    get_structure_graph_scene
from dash.dependencies import Output, Input
from dash_mp_components import Simple3DScene
from pydefect.analyzer.eigenvalue_plotter import EigenvaluePlotlyPlotter
from pydefect.corrections.efnv_correction.site_potential_plotter import \
    SitePotentialPlotlyPlotter

Legend.uniform_radius = 0.2


class SingleDefectComponent(MPComponent):

    default_title = "Defect structure"

    def __init__(self,
                 potential_plotter: SitePotentialPlotlyPlotter,
                 eigenvalue_plotter: EigenvaluePlotlyPlotter,
                 scene_dicts: SceneDicts,
                 name: str,
                 id: str = None,
                 scene_settings: Optional[Dict] = None,
                 **kwargs):

        self.scene_dicts = scene_dicts
        self.potential_plotter = potential_plotter
        self.eigenvalue_plotter = eigenvalue_plotter
        super().__init__(id=id, **kwargs)
        self.options = [{"label": i, "value": i} for i in scene_dicts.scenes.keys()]
        self.create_store("scene_dicts", initial_data=scene_dicts)

#        default_data=scene_dicts.structure_graph,
        self.initial_scene_settings = {"extractAxis": True,
                                       "renderer": "webgl",
                                       "defaultZoom": 3.0}
        if scene_settings:
            self.initial_scene_settings.update(scene_settings)

        # self.create_store("scene_settings",
        #                   initial_data=self.initial_scene_settings)

        initial_scene_additions = Scene(name="parchg", contents=[list(scene_dicts.scenes.values())[0]])
        self.create_store("scene_additions",
                          initial_data=initial_scene_additions)

        scene, legend = self.get_scene_and_legend(
            scene_dicts.structure_graph,
            scene_additions=initial_scene_additions)
        self.create_store("legend_data", initial_data=legend)
        self.create_store("graph", initial_data=scene_dicts.structure_graph)

        self.name = name

        # this is used by a Simple3DScene component, not a dcc.Store
        self._initial_data["scene"] = scene

    @property
    def _sub_layouts(self):
        pot_fig = self.potential_plotter.create_figure()
        eig_fig = self.eigenvalue_plotter.create_figure()
        pot_graph = dcc.Graph(figure=pot_fig)
        eig_graph = dcc.Graph(figure=eig_fig)

        choice_locpot = self.get_choice_input(kwarg_label="orbital",
                                              options=self.options,
                                              state={"orbital": self.options[0]["value"]})
        structure = html.Div(
            Simple3DScene(
                id=f"{self.id()}_scene",
                data=self.initial_data["scene"],
                settings=self.initial_scene_settings,
                sceneSize="100%",
            ),
            style={
                "width": "100%",
                "height": "100%",
                "overflow": "hidden",
                "margin": "0 auto",
            },
        )

        return {"pot_graph": pot_graph, "eig_graph": eig_graph,
                "locpot": choice_locpot, "structure": structure}

    def layout(self) -> html.Div:
        """
        :param size: a CSS string specifying width/height of Div
        :return: A html.Div containing the 3D structure or molecule
        """
        pot = html.Div(self._sub_layouts["pot_graph"])
        eig = html.Div(self._sub_layouts["eig_graph"])
        structure = html.Div([html.Div(self._sub_layouts["structure"], style={"width": "350px", "height": "350px"}),
                              self._sub_layouts["locpot"]])

        return Columns([H3(self.name), pot, eig, Column([structure])])

    @staticmethod
    def get_scene_and_legend(graph: StructureGraph, scene_additions=None
                             ) -> Tuple[Scene, Dict[str, str]]:

        legend = Legend(graph.structure,
                        color_scheme="VESTA",
                        radius_scheme="uniform",
                        cmap_range=None)

        scene = get_structure_graph_scene(
            graph,
            draw_image_atoms=True,
            bonded_sites_outside_unit_cell=False,
            hide_incomplete_edges=True,
            explicitly_calculate_polyhedra_hull=False,
            group_by_symmetry=False,
            legend=legend)

        scene.name = "DefectStructureComponentScene"

        # axes = graph.structure.lattice._axes_from_lattice()
        # axes.visible = True
        # scene.contents.append(axes)

        scene = scene.to_json()
        if scene_additions:
            # raise NotImplementedError
            scene["contents"].append(scene_additions)

        return scene, legend.get_legend()

    def generate_callbacks(self, app, cache):
        @app.callback(Output(self.id("scene"), "data"),
                      [Input(self.get_kwarg_id("orbital"), "value")])
        def update_structure(orbital):
            # orbital is ["0125_up"]
            scene_additions = Scene(name=orbital[0],
                                    contents=[self.scene_dicts.scenes[orbital[0]]])
                                              # Arrows(positionPairs=[[[0,0,0], [s,s,s]]],
                                              #        radius=0.1,
                                              #        headLength=0.6,
                                              #        clickable=True,
                                              #        headWidth=0.3)])
            scene, _ = self.get_scene_and_legend(
                self.scene_dicts.structure_graph,
                scene_additions=scene_additions)
            return scene
