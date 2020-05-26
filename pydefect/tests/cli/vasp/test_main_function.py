# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

import numpy as np
from monty.serialization import loadfn
from pymatgen import IStructure
from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.analyzer.unitcell import Unitcell
from pydefect.cli.vasp.main_function import make_supercell, make_defect_set, \
    make_defect_entries, make_unitcell
from pydefect.input_maker.defect import SimpleDefect
from pydefect.input_maker.defect_set import DefectSet


def test_make_unitcell(mocker):
    vasprun_band_mock = mocker.Mock(spec=Vasprun, autospec=True)
    outcar_band_mock = mocker.Mock(spec=Outcar, autospec=True)
    outcar_dielectric_mock = mocker.Mock(spec=Outcar, autospec=True)
    args = Namespace(vasprun_band=vasprun_band_mock,
                     outcar_band=outcar_band_mock,
                     outcar_dielectric=outcar_dielectric_mock)
    mock = mocker.patch("pydefect.cli.vasp.main_function.make_unitcell_from_vasp")
    mock.return_value = Unitcell(vbm=1.0,
                                 cbm=2.0,
                                 ele_dielectric_const=np.eye(3),
                                 ion_dielectric_const=np.eye(3))
    make_unitcell(args)
    mock.assert_called_once_with(vasprun_band=vasprun_band_mock,
                                 outcar_band=outcar_band_mock,
                                 outcar_dielectric=outcar_dielectric_mock)


def test_make_supercell_from_matrix(simple_cubic, simple_cubic_2x1x1, tmpdir):
    matrix = [2, 1, 1]
    args = Namespace(unitcell=simple_cubic, matrix=matrix, min_num_atoms=None, max_num_atoms=None)

    tmpdir.chdir()
    make_supercell(args)
    info = loadfn("supercell_info.json")
    assert IStructure.from_file("SPOSCAR") == simple_cubic_2x1x1
    assert info.structure == simple_cubic_2x1x1
    assert info.transform_matrix == [[2, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_make_recommended_supercell(simple_cubic, simple_cubic_2x2x2, tmpdir):
    args = Namespace(unitcell=simple_cubic, matrix=None, min_num_atoms=8, max_num_atoms=8)

    tmpdir.chdir()
    make_supercell(args)
    info = loadfn("supercell_info.json")
    assert IStructure.from_file("SPOSCAR") == simple_cubic_2x2x2
    assert info.structure == simple_cubic_2x2x2
    assert info.transform_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]


def test_make_defect_set(tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    args = Namespace(oxi_states={"He": 1}, dopants=["Li"], kwargs=["Li_H1", "Va_He1", "Va_H1_-1"])
    make_defect_set(args)

    simple_defects = {SimpleDefect(None, "He1", [-1, 0, 1]),
                      SimpleDefect(None, "H1", [-1]),
                      SimpleDefect("Li", "H1", [0]),
                      }
    assert DefectSet.from_yaml() == DefectSet(defects=simple_defects)


def test_make_defect_entries(tmpdir, supercell_info):
    tmpdir.chdir()
    supercell_info.to_json_file()
    defect_set = DefectSet({SimpleDefect(None, "He1", [-1, 0])})
    defect_set.to_yaml()
    args = Namespace()
    make_defect_entries(args)
    dir_names = {str(dirname) for dirname in Path(".").glob("*_*_*")}
    assert dir_names == {"Va_He1_-1", "Va_He1_0"}

    file_names = {str(file_name.name) for file_name in Path("Va_He1_-1").glob("*")}
    assert file_names == {"POSCAR", "defect_entry.json"}


def test_make_calc_results(tmpdir):
    args = Namespace(dirs=[Path("a"), Path("b")])
    Path("a").name


"""
TODO
- 

DONE
- make 2x2x2 supercell
- make supercell_info.json
- Generate defect_set from minimum args.
"""