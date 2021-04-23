# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pydefect.cli.vasp.make_unitcell import make_unitcell_from_vasp
from pymatgen.io.vasp import Vasprun, Outcar


def test_unitcell(vasp_files):
    """
    HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE, excluding Hartree and local field effects)
 ------------------------------------------------------
           1.269877     0.000000    -0.000000
           0.000000     1.269877     0.000000
           0.000000     0.000000     1.269877
 ------------------------------------------------------

 MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)
 ------------------------------------------------------
           1.255879     0.000000    -0.000000
          -0.000000     1.255879     0.000000
          -0.000000     0.000000     1.255879
 ------------------------------------------------------
    """

    path = vasp_files / "unitcell_Ne_solid"
    unitcell = make_unitcell_from_vasp(
        vasprun_band=Vasprun(path / "vasprun-band.xml"),
        outcar_band=Outcar(path / "OUTCAR-band"),
        outcar_dielectric_clamped=Outcar(path / "OUTCAR-dielectric"),
        outcar_dielectric_ionic=Outcar(path / "OUTCAR-dielectric"),
        )

    assert unitcell.system == "Ne"
    assert unitcell.vbm == -10.3168
    assert unitcell.cbm == 1.2042
    assert unitcell.ele_dielectric_const[0][0] == 1.255879
    assert unitcell.ion_dielectric_const[0][0] == 0.0



