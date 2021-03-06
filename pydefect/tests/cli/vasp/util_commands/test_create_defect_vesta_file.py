# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from pydefect.analyzer.defect_structure_comparator import SiteDiff
from pydefect.analyzer.defect_structure_info import DefectStructureInfo, \
    Displacement
from pydefect.cli.vasp.util_commands.create_defect_vesta_file import \
    make_defect_vesta_file
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

    def_str_info = DefectStructureInfo(initial_site_sym="3m",
                               final_site_sym="m",
                               site_diff=site_diff,
                               site_diff_from_initial=site_diff_from_init,
                               symprec=0.1,
                               dist_tol=0.2,
                               anchor_atom_idx=3,
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
    make_defect_vesta_file(defect_structure=structure, def_str_info=def_str_info, disp_enhance_factor=100)
