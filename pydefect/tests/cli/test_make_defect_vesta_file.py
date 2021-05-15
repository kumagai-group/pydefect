# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import DefectStructureInfo, \
    Displacement
from pydefect.cli.make_defect_vesta_file import \
    MakeDefectVestaFile
from pymatgen.core import Structure, Lattice


def test_defect_vesta_file(def_str_info, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    x = MakeDefectVestaFile(defect_structure_info=def_str_info, arrow_factor=100)
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
2 Li Li 1.0 0.500000 0.500000 0.500000
 0.0 0.0 0.0 
3 H H 1.0 0.250000 0.250000 0.250000
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
 1  Be1  0.5 94 215 123 94 215 123 204 0
 2  Li  0.5 134 224 116 134 224 116 204 0
 3  H  0.5 255 204 204 255 204 204 204 0
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
X3 Be2 H1

CELLP
10.000000 10.000000 10.000000 90.000000 90.000000 90.000000

STRUC
1 Be 2.1 1.0 0.500000 0.500000 0.500900
 0.0 0.0 0.0 
2 H H_i 1.0 0.270000 0.250000 0.250000
 0.0 0.0 0.0 
3 Be Be_Li 1.0 0.500000 0.500000 0.500900
 0.0 0.0 0.0 
4 XX Va_H 1.0 0.250000 0.250000 0.250000
 0.0 0.0 0.0 
5 XX Initial_H 1.0 0.250000 0.250000 0.250000
 0.0 0.0 0.0 
6 XX center 1.0 0.380000 0.375000 0.375175
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
 5  Initial_H  0.5 30 30 30 30 30 30 204 1
 6  center  0.3 30 30 30 30 30 30 204 0
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