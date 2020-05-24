# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.make_unitcell import make_unitcell_from_vasp


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

    path = vasp_files / "unitcell_He_solid"
    unitcell = make_unitcell_from_vasp(
        vasprun_band=Vasprun(path / "vasprun-band.xml"),
        outcar_band=Outcar(path / "OUTCAR-band"),
        outcar_dielectric=Outcar(path / "OUTCAR-dielectric"))

    assert unitcell.vbm == -10.3168
    assert unitcell.cbm == 1.2042
    assert unitcell.ele_dielectric_const[0][0] == 1.255879
    assert unitcell.ion_dielectric_const[0][0] == 0.0



"""
TODO
- Construct Unitcell only from dielectric calculations. 


    def set_band_edge_from_vasp(self,
                                directory_path: str,
                                vasprun_name: str = "vasprun.xml",
                                outcar_name: str = "OUTCAR") -> None:
        vasprun = Vasprun(os.path.join(directory_path, vasprun_name))
        outcar = Outcar(os.path.join(directory_path, outcar_name))

        # 2019/7/13 NEVER USE Vasprun.eigenvalue_band_properties
        # THERE IS A BUG TO ESTIMATE VBM AND CBM of lower band gap materials.
        _, vbm_info, cbm_info = band_gap_properties(vasprun, outcar)
        self.is_direct = vbm_info["kpoints"] == cbm_info["kpoints"]
        self._band_edge = [vbm_info["energy"], cbm_info["energy"]]
DONE
"""