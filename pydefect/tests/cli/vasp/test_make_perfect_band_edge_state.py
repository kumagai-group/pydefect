# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState, EdgeInfo, \
    OrbitalInfo
from pydefect.cli.vasp.make_perfect_band_edge_state import \
    make_perfect_band_edge_state_from_vasp
from pymatgen.io.vasp import Procar, Vasprun, Outcar


def test_make_perfect_band_edge_state_from_vasp(vasp_files):
    procar = Procar(vasp_files / "MgO_2x2x2_perfect" / "PROCAR")
    vasprun = Vasprun(vasp_files / "MgO_2x2x2_perfect" / "vasprun.xml", parse_potcar_file=False)
    outcar = Outcar(vasp_files / "MgO_2x2x2_perfect" / "OUTCAR")
    actual = make_perfect_band_edge_state_from_vasp(procar, vasprun, outcar)

    vbm_info = EdgeInfo(
        band_idx=127, kpt_coord=(0.25, 0.25, 0.25),
        orbital_info=OrbitalInfo(energy=2.7746,
                                 orbitals={'Mg': [0.0, 0.0, 0.0, 0.0],
                                           'O': [0.0, 0.704, 0.0, 0.0]},
                                 occupation=1.0))
    cbm_info = EdgeInfo(
        band_idx=128, kpt_coord=(0.25, 0.25, 0.25),
        orbital_info=OrbitalInfo(energy=8.2034,
                                 orbitals={'Mg': [0.192, 0.0, 0.0, 0.0],
                                           'O': [0.224, 0.096, 0.0, 0.0]},
                                 occupation=0.0))
    expected = PerfectBandEdgeState(vbm_info=vbm_info, cbm_info=cbm_info)
    assert actual == expected
