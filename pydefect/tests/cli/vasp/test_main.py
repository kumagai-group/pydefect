# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace

from pydefect.cli.vasp.main import parse_args


def test_make_supercell_wo_options(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.IStructure")
    parsed_args = parse_args(["ms"])
    # func is a pointer so need to point the same address.
    expected = Namespace(
        unitcell=mock.from_file.return_value,
        matrix=None,
        min_num_atoms=50,
        max_num_atoms=300,
        func=parsed_args.func,
    )
    assert parsed_args == expected
    mock.from_file.assert_called_once_with("POSCAR")


def test_make_supercell_w_options(mocker):
    mock = mocker.patch("pydefect.cli.vasp.main.IStructure")
    parsed_args = parse_args(["ms",
                              "-p", "POSCAR-tmp",
                              "--matrix", "1", "2", "3",
                              "--min_num_atoms", "1000",
                              "--max_num_atoms", "2000",
                              ]
                             )
    # func is a pointer so need to point the same address.
    expected = Namespace(
        unitcell=mock.from_file.return_value,
        matrix=[1, 2, 3],
        min_num_atoms=1000,
        max_num_atoms=2000,
        func=parsed_args.func,
    )
    assert parsed_args == expected
    mock.from_file.assert_called_once_with("POSCAR-tmp")


# def test_vasp_set_w_options():
#     parsed_args = parse_args(["vs",
#                               "--poscar", "POSCAR-tmp",
#                               "-t", "band",
#                               "-x", "pbesol",
#                               "-k", "4.2",
#                               "--potcar", "Mg_pv", "O_h",
#                               "-c", "10",
#                               "--user_incar_settings", "LREAD", "F",
#                               "-d", "c",
#                               "--vasprun", "vasprun_1",
#                               "--outcar", "OUTCAR_1",
#                               "--options", "encut", "800",
#                               "--uniform_kpt_mode",
#                               "--file_transfer_type", "WAVECAR", "C",
#                               ])
#
#     expected = Namespace(
#         poscar=Path("POSCAR-tmp"),
#         task=Task.band,
#         xc=Xc.pbesol,
#         kpt_density=4.2,
#         overridden_potcar=["Mg_pv", "O_h"],
#         charge=10.0,
#         user_incar_settings=["LREAD", "F"],
#         prev_dir=Path("c"),
#         vasprun=Path("vasprun_1"),
#         outcar=Path("OUTCAR_1"),
#         options=["encut", "800"],
#         uniform_kpt_mode=True,
#         file_transfer_type=["WAVECAR", "C"],
#         func=parsed_args.func,
#     )
#
#     assert parsed_args == expected
#

"""
TODO
-

DONE
"""