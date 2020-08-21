# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import Dict

import numpy as np
from crystal_toolkit.core.scene import Scene, Surface
from monty.json import MSONable
from pymatgen import Spin, Structure
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.io.vasp import Chgcar, VolumetricData
from skimage.measure import marching_cubes


@dataclass
class SceneDicts(MSONable):
    scene_dicts: Dict[str, Dict[str, np.ndarray]]
    structure_graph: StructureGraph

    @property
    def scenes(self):
        result = {}
        if len(self.scene_dicts) == 0:
            return {"empty_scene": Scene("empty", contents=[])}
        for k, v in self.scene_dicts.items():
            pos = [vert for triangle in v["vertices"][v["faces"]].tolist()
                   for vert in triangle]
            result[k] = Scene(k, contents=[Surface(positions=pos)])
        return result

    def as_dict(self) -> dict:
        """
        A JSON serializable dict representation of an object.
        """
        scene_dicts = {k: {k2: v2.tolist() for k2, v2 in v.items()}
                       for k, v in self.scene_dicts.items()}

        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "scene_dicts": scene_dicts,
                "structure_graph": self.structure_graph.as_dict()}

    @classmethod
    def from_dict(cls, d):
        scene_dicts = {k: {k2: np.array(v2) for k2, v2 in v.items()}
                       for k, v in d["scene_dicts"].items()}
        return cls(scene_dicts, StructureGraph.from_dict(d["structure_graph"]))


def make_scene_dicts(parchg_list, defect_pos):
    band_indices = [c.split(".")[-2] for c in parchg_list]
    # parchgs = [Chgcar.from_file(c) for c in parchg_list]
    parchgs = []
    for c in parchg_list:
        try:
            parchgs.append(Chgcar.from_file(c))
        except ValueError:
            continue
    if len(parchgs) == 0:
        return None
    structure: Structure = parchgs[0].structure
    if len(parchg_list) > 1:
        for parchg in parchgs[1:]:
            assert structure == parchg.structure

    lattice_matrix = parchgs[0].structure.lattice.matrix

    # centering to [0.5, 0.5, 0.5]
    shift_vector = np.array([0.5 - pos for pos in defect_pos])
    shape = parchgs[0].data["total"].shape
    shift_numbers = np.array(np.rint(shape * shift_vector), dtype=np.int)
    actual_shift_vector = shift_numbers / shape

    step_size = int(min(parchgs[0].dim) / 30)
    results = {}
    for name, parchg in zip(band_indices, parchgs):
        for spin in [Spin.up, Spin.down]:
            spin_str = "up" if spin is Spin.up else "down"
            key = name + "_" + spin_str
            data = parchg.spin_data[spin]
            x_shifted = np.roll(data, shift_numbers[0], axis=0)
            xy_shifted = np.roll(x_shifted, shift_numbers[1], axis=1)
            xyz_shifted = np.roll(xy_shifted, shift_numbers[2], axis=2)
            try:
                vertices, faces = get_vertices_and_faces(xyz_shifted,
                                                         lattice_matrix,
                                                         step_size=step_size)
                results[key] = {"vertices": vertices, "faces": faces}
            except RuntimeError:
                continue

    structure.translate_sites(indices=list(range(len(structure))),
                              vector=actual_shift_vector)

    graph = StructureGraph.with_empty_graph(structure=structure)
    return SceneDicts(results, structure_graph=graph)


def get_vertices_and_faces(array_data, lattice_matrix, step_size=10, level=None):
    # normals and values are removed.
    vertices, faces, _, _ = marching_cubes(array_data,
                                           level=level,
                                           step_size=step_size)
    # normalize to fractional coordinates
    vertices = vertices / array_data.shape
    # transform to cartesian
    vertices = np.dot(vertices - 0.5, lattice_matrix)
    return vertices, faces


@dataclass
class SingleScene(MSONable):
    scene_dict: Dict[str, np.ndarray]

    def scene(self, name=""):
        faces = self.scene_dict["faces"]
        pos = [vert for triangle in self.scene_dict["vertices"][faces].tolist()
               for vert in triangle]
        return Scene(name=name, contents=[Surface(positions=pos)])

    def as_dict(self) -> dict:
        scene_dict = {k: v.tolist() for k, v in self.scene_dict.items()}
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "scene_dict": scene_dict}

    @classmethod
    def from_dict(cls, d):
        scene_dicts = {k: np.array(v) for k, v in d["scene_dict"].items()}
        return cls(scene_dicts)


def make_single_scene(volmetric_data: VolumetricData, step_size=3):
    lattice_matrix = volmetric_data.structure.lattice.matrix
    data = volmetric_data.data["total"]
    vertices, faces = get_vertices_and_faces(data,
                                             lattice_matrix,
                                             step_size=step_size)
    return SingleScene({"vertices": vertices, "faces": faces})
