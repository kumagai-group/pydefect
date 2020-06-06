# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pydefect.analyzer.unitcell import Unitcell
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties


def make_unitcell_from_vasp(vasprun_band: Vasprun,
                            outcar_band: Outcar,
                            outcar_dielectric) -> Unitcell:
    outcar_dielectric.read_lepsilon()
    outcar_dielectric.read_lepsilon_ionic()
    band_edge_properties = VaspBandEdgeProperties(vasprun_band, outcar_band)
    vbm, cbm = band_edge_properties.vbm_cbm

    return Unitcell(
        vbm=vbm, cbm=cbm,
        ele_dielectric_const=outcar_dielectric.dielectric_tensor,
        ion_dielectric_const=outcar_dielectric.dielectric_ionic_tensor)
