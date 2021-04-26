# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path

from monty.serialization import loadfn
from pymatgen.core import Structure, Lattice

from pydefect.cli.vasp.make_poscars_from_query import make_poscars_from_query


def test_make_poscars_from_query(tmpdir):
    mg_structure = Structure(Lattice.cubic(1), species=["Mg", "Mg"],
                             coords=[[0, 0, 0], [0.5, 0.5, 0.5]])

    query = [{"full_formula": "Mg2", "task_id": "mp-1",
              "total_magnetization": 0.1,
              "band_gap": 0.0,
              "structure": mg_structure},
             {"full_formula": "O16"}]
    make_poscars_from_query(query, path=Path(tmpdir))
    tmpdir.chdir()
    s_mg = Structure.from_file(Path("Mg_mp-1") / "POSCAR")
    assert s_mg == mg_structure

    actual = loadfn(Path("Mg_mp-1") / "prior_info.yaml")
    expected = {'band_gap': 0.0, 'data_source': 'mp-1', 'total_magnetization': 0.1}
    assert actual == expected

    actual = Structure.from_file(Path("mol_O2") / "POSCAR")
    expected = Structure.from_str("""Default POSCAR of O2 molecule
   1.00000000000000     
    10.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   10.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   10.0000000000000000
   O 
     2
Direct
  0.0145089778962415  0.0000000000000000  0.0000000000000000
  0.1354910221037571  0.0000000000000000  0.0000000000000000""", fmt="POSCAR")
    assert actual == expected

