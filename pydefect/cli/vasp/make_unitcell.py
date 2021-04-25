# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pydefect.analyzer.unitcell import Unitcell
from pymatgen.io.vasp import Vasprun, Outcar
from vise.analyzer.vasp.band_edge_properties import VaspBandEdgeProperties


def make_unitcell_from_vasp(vasprun_band: Vasprun,
                            outcar_band: Outcar,
                            outcar_dielectric_clamped: Outcar,
                            outcar_dielectric_ionic: Outcar,
                            system_name: str = None) -> Unitcell:
    name = (system_name or
            vasprun_band.final_structure.composition.reduced_formula)
    outcar_dielectric_clamped.read_lepsilon()
    outcar_dielectric_ionic.read_lepsilon_ionic()
    band_edge_properties = VaspBandEdgeProperties(vasprun_band, outcar_band)
    vbm, cbm = band_edge_properties.vbm_cbm

    return Unitcell(
        system=name,
        # vbm and cbm are <class 'numpy.float64'>.
        vbm=float(vbm), cbm=float(cbm),
        ele_dielectric_const=outcar_dielectric_clamped.dielectric_tensor,
        ion_dielectric_const=outcar_dielectric_ionic.dielectric_ionic_tensor)
