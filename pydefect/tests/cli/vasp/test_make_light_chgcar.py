# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from pydefect.cli.vasp.make_light_chgcar import make_light_chgcar
from pymatgen import Structure, Lattice
from pymatgen.io.vasp import Chgcar
import numpy as np


def test_make_light_chgcar(tmpdir, mocker):
    print(tmpdir)
    tmpdir.chdir()
    structure = Structure(Lattice.cubic(1), species=["O"], coords=[[0]*3])
    chgcar = Chgcar(structure, data={"total": np.array([[[0.0, 0.1, -1e-5]]])})

    mock = mocker.patch("pydefect.cli.vasp.make_light_chgcar.add_density")

    make_light_chgcar(chgcar=chgcar, filename=Path("CHGCAR"),
                      vesta_file=Path("defect.vesta"),
                      to_vesta_file=Path("to.vesta"))

    actual = Chgcar.from_file("CHGCAR")
    expected = Chgcar(structure, data={"total": np.array([[[0, 9, 0]]])})
    np.testing.assert_almost_equal(actual.data["total"], expected.data["total"])

    # 0.13337 = 0.1 * 0.529 ** 3 * 9
    mock.assert_called_with(Path("defect.vesta"), Path("to.vesta"),
                            [0.13337, 0.66683, 1.06693], Path("CHGCAR"))
