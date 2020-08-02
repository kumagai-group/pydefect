# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.input_maker.add_interstitial import append_interstitial
from pydefect.input_maker.supercell_info import Interstitial
from pydefect.util.error_classes import NotPrimitiveError
from pymatgen import Structure, Lattice


def test_add_interstitial(cubic_supercell_info_wo_int):
    base_primitive = Structure.from_str("""
   1.00000000000000
    0 2.5 2.5
    2.5 0 2.5
    2.5 2.5 0
    H He
     1     1
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.5000000000000000  0.5000000000000000  0.5000000000000000""", fmt="POSCAR")
    new_supercell_info = append_interstitial(cubic_supercell_info_wo_int,
                                             base_primitive,
                                             [1/4, 1/4, 1/4])
    expected = Interstitial(frac_coords=[1/8, 1/8, 1/8], wyckoff_letter="d", site_symmetry="-43m")
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


"""
TODO
- check site symmetry


i1:
  representative_coords: [0.125, 0.125, 0.125]
  wyckoff: c
  site_symmetry: -43m
  multiplicity: 24
  coordination_distances:
    Mg: [1.82, 1.82, 1.82, 1.82]
    O: [1.82, 1.82, 1.82, 1.82]
  method: manual

DONE
- check primitive
- add interstitials to SupercellInfo
"""