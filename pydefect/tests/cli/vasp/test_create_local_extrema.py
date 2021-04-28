# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from pydefect.cli.vasp.create_local_extrema import \
    interstitials_from_volumetric_data, extrema_coords, \
    VolumetricDataAnalyzeParams, find_inequivalent_coords
from pydefect.input_maker.local_extrema import CoordInfo
from pydefect.util.structure_tools import Coordination
from pymatgen.io.vasp import Chgcar
import numpy as np


def test_interstitials_from_charge_density(vasp_files):
    aeccar0 = Chgcar.from_file(str(vasp_files / "NaMgF3_AECCAR0"))
    aeccar2 = Chgcar.from_file(str(vasp_files / "NaMgF3_AECCAR2"))
    aeccar = aeccar0 + aeccar2
    aeccar.write_file("CHGCAR")
    interstitials_from_volumetric_data(aeccar)


def test_find_inequivalent_coords(simple_cubic):
    actual = find_inequivalent_coords(simple_cubic,
                                      [np.array([0.1, 0.0, 0.0])],
                                      [1.0])
    expected = [CoordInfo((0.1, 0.0, 0.0), "4mm",
                          coordination=Coordination({"H": [0.1]}, 0.13, [0]),
                          quantity=1.0)]
    assert actual[0] == expected[0]


def test_extrema_coords(vasp_files):
    aeccar0 = Chgcar.from_file(str(vasp_files / "NaMgF3_AECCAR0"))
    aeccar2 = Chgcar.from_file(str(vasp_files / "NaMgF3_AECCAR2"))
    aeccar = aeccar0 + aeccar2
    actual = extrema_coords(volumetric_data=aeccar, find_min=True,
                            params=VolumetricDataAnalyzeParams())
    assert len(actual) == 36
    np.testing.assert_array_almost_equal(actual[0],
                                         np.array([0.2, 0.75, 0.908333]))

# def test_interstitials_from_potential():
#     locpot = Locpot.from_file(str(Path(__file__).parent / "vasp_files" / "NaMgF3_LOCPOT"))
#     interstitials_from_charge_density(locpot)

