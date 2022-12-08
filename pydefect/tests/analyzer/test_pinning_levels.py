# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pydefect.analyzer.defect_energy import DefectEnergySummary, DefectEnergies, \
    DefectEnergy
from pydefect.analyzer.pinning_levels import pinning_levels_from_charge_energies


def test_pinning_levels_from_charge_energies():
    defect_energies = DefectEnergies(
            atom_io={"O": -1},
            charges=[0, 1, 2],
            defect_energies=[DefectEnergy(1.0, {"corr": 2.0}, is_shallow=False),
                             DefectEnergy(2.0, {"corr": 2.0}, is_shallow=False),
                             DefectEnergy(3.0, {"corr": 2.0}, is_shallow=True)])

    defect_energy_summary = DefectEnergySummary(
        title="test",
        defect_energies={"Va_O1": defect_energies},
        rel_chem_pots={"A": {"O": -5.0}},
        cbm=2.0, supercell_vbm=-1.0, supercell_cbm=3.0)
    cbm = defect_energy_summary.cbm
    charge_energies = defect_energy_summary.charge_energies("A",
                                                            allow_shallow=True,
                                                            with_corrections=True,
                                                            e_range=(0, cbm))
    actual = pinning_levels_from_charge_energies(charge_energies=charge_energies)
    expected = """defect    hole pinning          electron pinning
--------  --------------------  ------------------
Va_O      charge 1, level 1.00  -"""
    assert actual == expected