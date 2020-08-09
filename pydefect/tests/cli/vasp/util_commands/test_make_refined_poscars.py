# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from pydefect.cli.vasp.util_commands.make_refined_poscar import \
    make_refined_structure
from pydefect.defaults import defaults
from pydefect.input_maker.defect_entry import DefectEntry
from pymatgen import Structure


def test_make_refined_structure(tmpdir, mocker, simple_cubic):
    tmpdir.chdir()
    tmpdir.join("Va_O1_2").mkdir()
    simple_cubic.to(fmt="POSCAR", filename=defaults.contcar)
    tmpdir.join(defaults.outcar).write("")
    tmpdir.join(defaults.vasprun).write("")

    mock_defect_entry = mocker.Mock(spec=DefectEntry, autospec=True)

    def side_effect(key):
        if str(key) == "defect_entry.json":
            mock_defect_entry.anchor_atom_index = 0
            mock_defect_entry.structure = simple_cubic
            return mock_defect_entry
        else:
            raise ValueError

    mock_loadfn = mocker.patch("pydefect.cli.vasp.util_commands.make_refined_poscar.loadfn", side_effect=side_effect)
    mock = mocker.patch("pydefect.cli.vasp.util_commands.make_refined_poscar.symmetrize_defect_structure")
    mock.return_value = simple_cubic

    make_refined_structure()

    mock_loadfn.assert_any_call("defect_entry.json")

    assert Structure.from_file(Path("CONTCAR.sym_1")) == simple_cubic
    assert Structure.from_file(Path("POSCAR")) == simple_cubic
