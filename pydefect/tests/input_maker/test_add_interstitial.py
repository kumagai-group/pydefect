# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.input_maker.add_interstitial import append_interstitial
from pydefect.input_maker.supercell_info import Interstitial
from pydefect.util.error_classes import NotPrimitiveError
from pymatgen import Structure, Lattice


def test_add_interstitial(cubic_supercell_info_wo_int):
    primitive = Structure(Lattice.rhombohedral(7.071068, 60),
                          species=["H", "He"],
                          coords=[[0.0]*3, [0.5]*3])
    new_supercell_info = append_interstitial(cubic_supercell_info_wo_int,
                                             primitive,
                                             [1/4, 1/4, 1/4])
    expected = Interstitial(frac_coords=[1/8, 1/8, 1/8],
                            wyckoff_letter="d",
                            site_symmetry="-43m")
    assert new_supercell_info.interstitials[0] == expected


def test_add_interstitial_not_primitive_error(cubic_supercell_info_wo_int):
    conventional_cell = Structure(Lattice.cubic(10),
                                  species=["H"] * 4 + ["He"] * 4,
                                  coords=[[0.0, 0.0, 0.0],
                                          [0.5, 0.5, 0.0],
                                          [0.5, 0.0, 0.5],
                                          [0.0, 0.5, 0.5],

                                          [0.0, 0.0, 0.5],
                                          [0.0, 0.5, 0.0],
                                          [0.5, 0.0, 0.0],
                                          [0.5, 0.5, 0.5],
                                        ])
    with pytest.raises(NotPrimitiveError):
        append_interstitial(cubic_supercell_info_wo_int, conventional_cell,
                            [1/4, 1/4, 1/4])

