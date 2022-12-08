# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from pydefect.analyzer.defect_energy import ChargeEnergies
from tabulate import tabulate


def pinning_levels_from_charge_energies(charge_energies: ChargeEnergies):

    result = []
    for name, charge_energy in charge_energies.charge_energies_dict.items():
        pin_levels = charge_energy.pinning_level(0, charge_energies.e_max)
        h_pin, e_pin = "-", "-"
        if pin_levels[0]:
            h_pin = f"charge {pin_levels[0][1]}, level {pin_levels[0][0]:.2f}"
        if pin_levels[1]:
            e_pin = f"charge {pin_levels[1][1]}, level {pin_levels[1][0]:.2f}"
        result.append([name, h_pin, e_pin])

    headers = ["defect", "hole pinning", "electron pinning"]
    return tabulate(result, headers=headers)