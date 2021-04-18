# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

from pydefect.analyzer.defect_energy import DefectEnergyInfo, DefectEnergy, \
    DefectEnergySummary, DefectEnergies
from pydefect.analyzer.make_defect_energy_summary import \
    make_defect_energy_summary


def test_make_defect_energy_summary(mocker):
    de1 = DefectEnergy(0.0, {"PC correction": 2.0}, False)
    de2 = DefectEnergy(1.0, {"PC correction": 3.0}, True)
    defect_infos = [DefectEnergyInfo("Va_Mg1", 0, {"Mg": -1}, de1),
                    DefectEnergyInfo("Va_Mg1", 1, {"Mg": -1}, de2)]
    rel_chem_pot = {"A": {"Mg": 5.0}}

    unitcell = mocker.Mock()
    perf_be_state = mocker.Mock()

    actual = make_defect_energy_summary(defect_infos, rel_chem_pot, unitcell,
                                        perf_be_state)

    defect_energies = {"Va_Mg1": DefectEnergies(atom_io={"Mg": -1},
                                                charges=[0, 1],
                                                defect_energies=[de1, de2])}
    expected = DefectEnergySummary(title=unitcell.system,
                                   defect_energies=defect_energies,
                                   rel_chem_pots=rel_chem_pot,
                                   cbm=unitcell.cbm,
                                   supercell_vbm=perf_be_state.vbm_info.energy,
                                   supercell_cbm=perf_be_state.cbm_info.energy)
    assert actual == expected
