# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

import fire
from monty.serialization import loadfn
from pymatgen import Structure, Lattice
from pymatgen.io.vasp import Kpoints
from pymatgen.io.vasp.sets import MPRelaxSet
from vise.input_set.datasets.potcar_set import PotcarSet
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.task import Task
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.input_set.xc import Xc

nupdown = loadfn(Path(__file__).parent / "atom_magnetization.yaml")


def make_atom_vasp_set(potcar_set: PotcarSet, xc: Xc):
    for element in PotcarSet.mp_relax_set.potcar_dict().keys():
        Path(element).mkdir()
        structure = Structure(Lattice.cubic(10),
                              coords=[[0.5]*3], species=[element])
        input_options = CategorizedInputOptions(structure,
                                                task=Task.cluster_opt,
                                                xc=xc,
                                                potcar_set=potcar_set)
        vasp_input_files = VaspInputFiles(
            input_options,
            overridden_incar_settings={"ISPIN": 2,
                                       "NUPDOWN": nupdown[element],
                                       "NELM": 300})
        vasp_input_files.create_input_files(dirname=Path(element))


def make_atom_mp_relax_set():
    for element, potcar in PotcarSet.mp_relax_set.potcar_dict().items():
        if potcar is None:
            continue
        Path(element).mkdir()
        structure = Structure(Lattice.cubic(10),
                              coords=[[0.5]*3], species=[element])

        mp_set = MPRelaxSet(structure,
                            user_kpoints_settings=Kpoints(kpts=((1, 1, 1),)),
                            user_incar_settings={"ISIF": 2,
                                                 "ISMEAR": 0,
                                                 "NUPDOWN": nupdown[element],
                                                 "NELM": 300})
        mp_set.write_input(element)


if __name__ == '__main__':
    def make_set(_set, xc):
        make_atom_vasp_set(PotcarSet.from_string(_set), Xc.from_string(xc))

    fire.Fire(make_set)
