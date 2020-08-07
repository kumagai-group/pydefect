# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import shutil

import fire
from monty.serialization import loadfn
from pydefect.analyzer.defect_structure_analyzer import \
    symmetrize_defect_structure
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Structure
from vise.util.structure_symmetrizer import StructureSymmetrizer


def make_refined_structure():
    defect_entry: DefectEntry = loadfn("defect_entry.json")
    structure = Structure.from_file(defaults.contcar)
    symmetrizer = StructureSymmetrizer(structure,
                                       defaults.symmetry_length_tolerance,
                                       defaults.symmetry_angle_tolerance)
    if symmetrizer.point_group == "1":
        return

    shutil.move(defaults.contcar, str(defaults.contcar) + ".sym_1")
    shutil.move(defaults.outcar, str(defaults.outcar) + ".sym_1")
    shutil.move(defaults.vasprun, str(defaults.vasprun) + ".sym_1")

    refined_structure = symmetrize_defect_structure(
        symmetrizer,
        defect_entry.anchor_atom_index,
        defect_entry.anchor_atom_coords)

    refined_structure.to(fmt="POSCAR", filename="POSCAR")


if __name__ == '__main__':
    fire.Fire(make_refined_structure)
