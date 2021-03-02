# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

import numpy as np
import pytest
from crystal_toolkit.core.scene import Surface, Scene
from pydefect.analyzer.dash_components.scenes_from_volumetric_data \
    import SceneDicts, get_vertices_and_faces, make_scene_dicts, SingleScene
from pymatgen import Structure
from pymatgen.analysis.local_env import MinimumDistanceNN
from vise.tests.helpers.assertion import assert_msonable


@pytest.fixture
def scene_dicts():
    structure = Structure.from_file(str(Path(__file__).parent / "mgo_defects" / "Va_O1_0" / "CONTCAR-finish"))
    graph = MinimumDistanceNN().get_bonded_structure(structure=structure)
    return SceneDicts(
        {"125_up": {"vertices": np.array([[0] * 3, [0.1] * 3, [0.2] * 3]),
                    "faces": np.array([[0, 1, 2]])}},
        structure_graph=graph)


@pytest.fixture
def scene_dict():
    return SingleScene(
        {"vertices": np.array([[0] * 3, [0.1] * 3, [0.2] * 3]),
         "faces": np.array([[0, 1, 2]])})


def test_scene_dicts_msonable(scene_dicts, scene_dict):
    assert_msonable(scene_dicts)
    assert_msonable(scene_dict)


def test_scenes(scene_dicts, scene_dict):
    expected = Scene("125_up", contents=[Surface(positions=[[0]*3, [0.1]*3, [0.2]*3])])
    assert scene_dicts.scenes == {"125_up": expected}
    assert scene_dict.scene("125_up") == expected


def test_make_scene_dicts(mocker):
    actual = make_scene_dicts([str(Path(__file__).parent / "mgo_defects" / "actual_parchg_Va_O1_2" / "PARCHG.0125.ALLK")], [0.25, 0, 0])
    assert list(actual.scene_dicts.keys()) == ["0125_up", "0125_down"]

    mock = mocker.patch("pydefect.analyzer.dash_components.scenes_from_volumetric_data.Chgcar")
    mock.from_file.side_effect = ValueError
    actual = make_scene_dicts([str(Path(__file__).parent / "mgo_defects" / "actual_parchg_Va_O1_2" / "PARCHG.0125.ALLK")], [0.25, 0, 0])
    assert actual is None


def test_get_vertices_and_faces(mocker):
    mock = mocker.patch("pydefect.analyzer.dash_components.scenes_from_volumetric_data.marching_cubes", autospec=True)
    vertices = np.array([[0] * 3, [0.1] * 3, [0.2] * 3], dtype=np.float32)
    faces = np.array([[0, 1, 2]], dtype=np.int32)
    mock.return_value = vertices, faces, None, None
    array_data = np.array([[[0.1, 0.2, 0.3, 0.4], [0.5, 0.6, 0.7, 0.8]]])  # shape is (1, 2, 4)
    lattice_matrix = np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]], dtype=np.float32)
    actual_vertices, actual_faces = get_vertices_and_faces(array_data, lattice_matrix, level=10, step_size=20)
    expected_vertices = np.array([[-5.0, -5.0, -5.0],
                                  [-4.0, -4.5, -4.75],
                                  [-3.0, -4.0, -4.5]])
    mock.assert_called_once_with(array_data, step_size=20, level=10 * 0.8)

    np.testing.assert_array_almost_equal(actual_vertices, expected_vertices)
    np.testing.assert_array_almost_equal(actual_faces, faces)


