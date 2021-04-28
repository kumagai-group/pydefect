# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import pytest
from pydefect.input_maker.local_extrema import VolumetricDataLocalExtrema, \
    CoordInfo, VolumetricDataAnalyzeParams
from pydefect.util.structure_tools import Coordination
from vise.tests.helpers.assertion import assert_json_roundtrip


@pytest.fixture
def local_extrema(simple_cubic):
    coordination = Coordination({"Mn": [1.0, 2.0]}, cutoff=4.0,
                                neighboring_atom_indices=[1, 2])
    local_extremum = CoordInfo(site_symmetry="1",
                               coordination=coordination,
                               frac_coords=[(0.1, 0.1, 0.1)],
                               quantities=[2.1])

    return VolumetricDataLocalExtrema(unit_cell=simple_cubic, is_min=True,
                                      extrema_points=[local_extremum],
                                      info="test",
                                      params=VolumetricDataAnalyzeParams())


def test_local_extrema_json_roundtrip(local_extrema, tmpdir):
    tmpdir.chdir()
    assert_json_roundtrip(local_extrema, tmpdir)


