# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from itertools import groupby
from typing import List, Dict

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from pydefect.analyzer.defect_energy import DefectEnergyInfo, \
    DefectEnergySummary, DefectEnergies
from pydefect.analyzer.unitcell import Unitcell
from pymatgen import Element


def make_defect_energy_summary(
        single_energies: List[DefectEnergyInfo],
        rel_chem_pot: Dict[str, Dict[Element, float]],
        unitcell: Unitcell,
        perf_be_state: PerfectBandEdgeState) -> DefectEnergySummary:

    defect_energies = {}
    for _, grouped_es in groupby(single_energies, lambda x: x.name):
        grouped_es = list(grouped_es)
        name, atom_io = grouped_es[0].name, grouped_es[0].atom_io
        charges, des = [], []
        for es in grouped_es:
            assert atom_io == es.atom_io
            charges.append(es.charge)
            des.append(es.defect_energy)
        defect_energies[name] = DefectEnergies(atom_io, charges, des)

    return DefectEnergySummary(title=unitcell.title,
                               defect_energies=defect_energies,
                               rel_chem_pots=rel_chem_pot,
                               cbm=unitcell.cbm,
                               supercell_vbm=perf_be_state.vbm_info.energy,
                               supercell_cbm=perf_be_state.cbm_info.energy)
