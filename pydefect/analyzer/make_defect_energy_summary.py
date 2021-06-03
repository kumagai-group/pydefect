# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from itertools import groupby
from typing import List, Dict

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from pydefect.analyzer.defect_energy import DefectEnergyInfo, \
    DefectEnergySummary, DefectEnergies
from pydefect.analyzer.unitcell import Unitcell
from pydefect.chem_pot_diag.chem_pot_diag import TargetVertices


def make_defect_energy_summary(
        energy_infos: List[DefectEnergyInfo],
        target_vertices: TargetVertices,
        unitcell: Unitcell,
        perfect_band_edge: PerfectBandEdgeState) -> DefectEnergySummary:

    defect_energies = {}
    key = lambda x: x.name
    # MUST NEED SORTED.
    for _, grouped_es in groupby(sorted(energy_infos, key=key), key=key):
        grouped_es = list(grouped_es)
        name, atom_io = grouped_es[0].name, grouped_es[0].atom_io
        charges, des = [], []
        for es in grouped_es:
            assert atom_io == es.atom_io
            charges.append(es.charge)
            des.append(es.defect_energy)
        defect_energies[name] = DefectEnergies(atom_io, charges, des)

    return DefectEnergySummary(title=unitcell.system,
                               defect_energies=defect_energies,
                               rel_chem_pots=target_vertices.chem_pots,
                               cbm=unitcell.cbm - unitcell.vbm,
                               supercell_vbm=perfect_band_edge.vbm_info.energy - unitcell.vbm,
                               supercell_cbm=perfect_band_edge.cbm_info.energy - unitcell.vbm)
