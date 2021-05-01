# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from pydefect.analyzer.defect_energy import DefectEnergyInfo, DefectEnergy, \
    DefectEnergySummary, DefectEnergies
from pydefect.analyzer.make_defect_energy_summary import \
    make_defect_energy_summary
from pydefect.chem_pot_diag.chem_pot_diag import TargetVertices, TargetVertex


def test_make_defect_energy_summary(mocker):
    energy1 = DefectEnergy(0.0, {"PC correction": 2.0}, False)
    energy2 = DefectEnergy(1.0, {"PC correction": 3.0}, True)
    defect_infos = [DefectEnergyInfo("Va_Mg1", 0, {"Mg": -1}, energy1),
                    DefectEnergyInfo("Va_Mg1", 1, {"Mg": -1}, energy2)]
    target_vertices = TargetVertices(
        target="MgO", vertices={"A": TargetVertex({"Mg": 5.0})})

    unitcell = mocker.Mock()
    unitcell.vbm = 1.0
    unitcell.cbm = 11.0

    perf_be_state = mocker.Mock()
    perf_be_state.vbm_info.energy = 0.0
    perf_be_state.cbm_info.energy = 12.0

    actual = make_defect_energy_summary(defect_infos, target_vertices, unitcell,
                                        perf_be_state)

    defect_energies = {"Va_Mg1": DefectEnergies(
        atom_io={"Mg": -1}, charges=[0, 1], defect_energies=[energy1, energy2])}

    expected = DefectEnergySummary(title=unitcell.system,
                                   defect_energies=defect_energies,
                                   rel_chem_pots={"A": {"Mg": 5.0}},
                                   cbm=10.0,
                                   supercell_vbm=-1.0,
                                   supercell_cbm=11.0)
    assert actual == expected
