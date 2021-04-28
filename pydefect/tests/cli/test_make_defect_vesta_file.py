# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import DefectStructureInfo, \
    Displacement
from pydefect.cli.make_defect_vesta_file import \
    MakeDefectVestaFile
from pymatgen.core import Structure, Lattice


def test_defect_vesta_file(tmpdir):
    print(tmpdir)
    site_diff = SiteDiff(removed=[(0, 'H', (0.25, 0.25, 0.25))],
                         inserted=[(0, 'H', (0.27, 0.25, 0.25))],
                         removed_by_sub=[(2, 'Li', (0.5, 0.5, 0.5))],
                         inserted_by_sub=[(2, 'Be', (0.5, 0.5, 0.5009))])
    site_diff_from_init = SiteDiff(removed=[(0, 'H', (0.25, 0.25, 0.25))],
                                   inserted=[(0, 'H', (0.27, 0.25, 0.25))],
                                   removed_by_sub=[(2, 'Li', (0.5, 0.5, 0.5))],
                                   inserted_by_sub=[(2, 'Be', (0.5, 0.5, 0.5))])
    displacements = [None,
                     Displacement(specie="He",
                                  original_pos=(0.75, 0.75, 0.75),
                                  final_pos=(0.76, 0.75, 0.75),
                                  distance_from_defect=6.465437190360757,
                                  disp_vector=(0.1, 0.0, 0.0),
                                  displace_distance=0.09999999999999964,
                                  angle=124.9),
                     Displacement(specie="Be",
                                  original_pos=(0.5, 0.5, 0.5),
                                  final_pos=(0.5, 0.5, 0.5009),
                                  distance_from_defect=2.1355627039494767,
                                  disp_vector=(0.0, 0.0, 0.009),
                                  displace_distance=0.009000000000000341,
                                  angle=125.8),
                     Displacement(specie="U",
                                  original_pos=(0, 0, 0),
                                  final_pos=(0.0, 0.0, 0.0),
                                  distance_from_defect=6.525191802736529,
                                  disp_vector=(0.0, 0.0, 0.0),
                                  displace_distance=0.0,
                                  angle=None)]

    def_str_info = DefectStructureInfo(
        initial_site_sym="3m",
        final_site_sym="m",
        site_diff=site_diff,
        site_diff_from_initial=site_diff_from_init,
        symprec=0.1,
        dist_tol=0.2,
        anchor_atom_idx=3,
        neighbor_atom_indices=[1],
        neighbor_cutoff_factor=1.0000,
        drift_vector=(0.0, 0.0, 0.0001),
        drift_dist=0.001,
        center=(0.38, 0.375, 0.37517500000000004),
        displacements=displacements)
    structure = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.27, 0.25, 0.2501],
                [0.76, 0.75, 0.7501],
                [0.5, 0.5, 0.501],
                [0.0, 0.0, 0.0001]])
    tmpdir.chdir()
    x = MakeDefectVestaFile(defect_structure=structure, defect_structure_info=def_str_info, arrow_factor=100)
    x.initial_vesta.write_file("defect_initial.vesta")
    x.final_vesta.write_file("defect.vesta")

    expected_initial = """#VESTA_FORMAT_VERSION 3.5.0
TITLE
Li1 Be1 H1

CELLP
10.000000 10.000000 10.000000 90.000000 90.000000 90.000000

STRUC
1 Be Be1 1.0 0.500000 0.500000 0.500000
 0.0 0.0 0.0 
2 Li Li2 1.0 0.500000 0.500000 0.500000
 0.0 0.0 0.0 
3 H H3 1.0 0.250000 0.250000 0.250000
 0.0 0.0 0.0 
 0 0 0 0 0 

BOUND
-0.120000 0.880000 -0.125000 0.875000 -0.124825 0.875175
 0 0 0 0 0 

SBOND
1 Li Be 0.0 1.788  0  0  1  0  1 0.12   2.000 161  33 246
2 Li H 0.0  1.08  0  0  1  0  1 0.12   2.000 161  33 246
3 Be Li 0.0 1.788  0  0  1  0  1 0.12   2.000 161  33 246
4 Be H 0.0 0.708  0  0  1  0  1 0.12   2.000 161  33 246
5 H Li 0.0  1.08  0  0  1  0  1 0.12   2.000 161  33 246
6 H Be 0.0 0.708  0  0  1  0  1 0.12   2.000 161  33 246
 0 0 0 0 

SITET
 1  Be1  0.5 94 215 123 94 215 123 204 1
 2  Li2  0.5 134 224 116 134 224 116 204 1
 3  H3  0.5 255 204 204 255 204 204 204 1
 0 0 0 0 



STYLE
VECTS  1.0
SECTS  64  1
SECTP 
   1 0 0  0.00000E+00  0.00000E+00  -1.00000E+01  1.00000E+01
UCOLP 
 0  2  1.000   0   0   0
ATOMS  1  0  1
BONDP 
   1  16  0.12  1.000 127 127 127"""
    assert Path("defect_initial.vesta").read_text() == expected_initial

    expected_final = """#VESTA_FORMAT_VERSION 3.5.0
TITLE
X2 Be2 H1

CELLP
10.000000 10.000000 10.000000 90.000000 90.000000 90.000000

STRUC
1 Be 2.1 1.0 0.500000 0.500000 0.501000
 0.0 0.0 0.0 
2 H H_i 1.0 0.270000 0.250000 0.250000
 0.0 0.0 0.0 
3 Be Be_Li 1.0 0.500000 0.500000 0.500900
 0.0 0.0 0.0 
4 XX Va_H 1.0 0.250000 0.250000 0.250000
 0.0 0.0 0.0 
5 XX center 1.0 0.380000 0.375000 0.375175
 0.0 0.0 0.0 
 0 0 0 0 0 

BOUND
-0.120000 0.880000 -0.125000 0.875000 -0.124825 0.875175
 0 0 0 0 0 

SBOND
4 Be H 0.0 0.708  0  0  1  0  1 0.12   2.000 161  33 246
6 H Be 0.0 0.708  0  0  1  0  1 0.12   2.000 161  33 246
 0 0 0 0 

SITET
 1  2.1  0.5 94 215 123 94 215 123 204 1
 2  H_i  0.5 255 204 204 255 204 204 204 1
 3  Be_Li  0.5 94 215 123 94 215 123 204 1
 4  Va_H  0.5 30 30 30 30 30 30 204 1
 5  center  0.3 30 30 30 30 30 30 204 0
 0 0 0 0 



STYLE
VECTS  1.0
SECTS  64  1
SECTP 
   1 0 0  0.00000E+00  0.00000E+00  -1.00000E+01  1.00000E+01
UCOLP 
 0  2  1.000   0   0   0
ATOMS  1  0  1
BONDP 
   1  16  0.12  1.000 127 127 127"""

    assert Path("defect.vesta").read_text() == expected_final