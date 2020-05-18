# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from pydefect.input_maker.defect_entries_maker import DefectEntriesMaker
from pydefect.input_maker.defect_set import DefectSet
from pydefect.input_maker.defect_set_maker import DefectSetMaker
from pydefect.input_maker.supercell_info import SupercellInfo
from pydefect.input_maker.supercell_maker import SupercellMaker


def make_supercell(args):
    if args.matrix:
        maker = SupercellMaker(args.unitcell, args.matrix)
    else:
        kwargs = {}
        if args.min_num_atoms:
            kwargs["min_num_atoms"] = args.min_num_atoms
        if args.max_num_atoms:
            kwargs["max_num_atoms"] = args.max_num_atoms
        maker = SupercellMaker(args.unitcell, **kwargs)

    maker.supercell.structure.to(filename="SPOSCAR")
    maker.supercell_info.to_json_file()


def make_defect_set(args):
    supercell_info = SupercellInfo.from_json_file()
    maker = DefectSetMaker(supercell_info,
                           args.oxi_states,
                           args.dopants,
                           keywords=args.kwargs)
    maker.defect_set.to_yaml()


def make_defect_entries(args):
    supercell_info = SupercellInfo.from_json_file()
    defect_set = DefectSet.from_yaml()
    maker = DefectEntriesMaker(supercell_info, defect_set)
    for defect_entry in maker.defect_entries:
        dir_path = Path(defect_entry.full_name)
        dir_path.mkdir()
        defect_entry.perturbed_structure.to(filename=dir_path / "POSCAR")
        defect_entry.to_json_file(filename=dir_path / "defect_entry.json")
