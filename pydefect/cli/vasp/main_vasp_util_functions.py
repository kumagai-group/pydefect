# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pydefect.cli.vasp.make_composition_energies_from_mp import \
    make_composition_energies_from_mp


def composition_energies_from_mp(args) -> None:
    composition_energies = make_composition_energies_from_mp(
        elements=args.elements, atom_energy_yaml=args.atom_energy_yaml)
    composition_energies.to_yaml_file()

