# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import fire
from pydefect.input_maker.defect_entry import make_defect_entry
from pymatgen import IStructure
from vise.util.logger import get_logger

logger = get_logger(__name__)


def make_defect_entry_from_poscars(dirname, perf_poscar, defect_poscar):
    names = dirname.split("_")
    try:
        name = "_".join(names[:-1])
        charge = int(names[-1])
    except ValueError:
        logger.warning(f"Directory name {dirname} is invalid.")
        raise

    perfect = IStructure.from_file(perf_poscar)
    defect = IStructure.from_file(defect_poscar)

    defect_entry = make_defect_entry(name=name,
                                     charge=charge,
                                     perfect_structure=perfect,
                                     defect_structure=defect)
    defect_entry.to_json_file()


if __name__ == '__main__':
    fire.Fire(make_defect_entry_from_poscars)
