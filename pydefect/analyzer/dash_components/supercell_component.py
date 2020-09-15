# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from typing import Tuple

from crystal_toolkit.core.legend import Legend
from crystal_toolkit.core.mpcomponent import MPComponent
from crystal_toolkit.core.scene import Scene
from crystal_toolkit.helpers.layouts import *
from crystal_toolkit.renderables.site import Site
from crystal_toolkit.renderables.structuregraph import get_structure_graph_scene
from dash_mp_components import Simple3DScene
from pydefect.input_maker.supercell_info import SupercellInfo
from pymatgen import Structure, DummySpecie, Lattice
from pymatgen.analysis.graphs import StructureGraph


class SupercellComponent(MPComponent):
    def __init__(self, supercell_info: SupercellInfo, *args, **kwargs):
        s = Structure.from_dict(supercell_info.structure.as_dict())
        coords = [[f % 1 for f in x.frac_coords] for x in s]
        structure = Structure(lattice=s.lattice, species=s.species,
                              coords=coords)
        self.graph = StructureGraph.with_empty_graph(structure=structure)

        self.transformation_matrix = supercell_info.transformation_matrix
        self.interstitials = supercell_info.interstitials
        super().__init__(*args, **kwargs)

        self.title = "supercell size: " + supercell_info._transform_matrix_str

        scene, legend = self.get_scene_and_legend()
        self.create_store("legend_data", initial_data=legend)
        self.create_store("graph", initial_data=self.graph)

        # this is used by a Simple3DScene component, not a dcc.Store
        self._initial_data["scene"] = scene

    def get_scene_and_legend(self) -> Tuple[Scene, Dict[str, str]]:

        legend = Legend(self.graph.structure,
                        color_scheme="VESTA",
                        cmap_range=None)
        legend.uniform_radius = 0.3

        scene = get_structure_graph_scene(
            self.graph,
            draw_image_atoms=True,
            bonded_sites_outside_unit_cell=False,
            hide_incomplete_edges=True,
            explicitly_calculate_polyhedra_hull=False,
            group_by_symmetry=True,
            draw_polyhedra=False,
            legend=legend)

        scene.name = "DefectStructureComponentScene"

        lattice: Lattice = self.graph.structure.lattice
        origin = -self.graph.structure.lattice.get_cartesian_coords([0.5, 0.5, 0.5])
        scene_json = scene.to_json()
        for idx, i in enumerate(self.interstitials, 1):
            site = Site(species=DummySpecie(),
                        coords=lattice.get_cartesian_coords(i.frac_coords))
            interstitial_scene = site.get_scene(origin=origin)
            interstitial_scene.name = f"i{idx}"
            interstitial_scene.contents[0].contents[0].tooltip = f"i{idx}"
            scene_json["contents"].append(interstitial_scene.to_json())

        return scene_json, legend.get_legend()

    @property
    def _sub_layouts(self):
        structure = html.Div(
            Simple3DScene(
                id=f"{self.id()}_scene",
                data=self.initial_data["scene"],
                sceneSize="100%",
            ),
            style={
                "width": "600px",
                "height": "600px",
                "overflow": "hidden",
                "margin": "0 auto",
            },
        )

        return {"title": H4(self.title), "structure": structure}

    @property
    def layout(self) -> html.Div:
        """
        :param size: a CSS string specifying width/height of Div
        :return: A html.Div containing the 3D structure or molecule
        """
        return html.Div(Column([self._sub_layouts["title"],
                                self._sub_layouts["structure"]]))

    def generate_callbacks(self, app, cache):
        pass

