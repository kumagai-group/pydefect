# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    CompositionEnergy
from pydefect.cli.vasp.main_vasp_util_functions import \
    composition_energies_from_mp
from pymatgen.core import Composition


def test_composition_energies_from_mp(mocker, tmpdir):
    tmpdir.chdir()
    args = Namespace(elements=["H"], atom_energy_yaml="a.yaml")
    mock = mocker.patch("pydefect.cli.vasp.main_vasp_util_functions.make_composition_energies_from_mp")
    mock.return_value = CompositionEnergies(
        {Composition("H"): CompositionEnergy(1.0, "mp-1")})
    composition_energies_from_mp(args)
    mock.assert_called_once_with(elements=["H"], atom_energy_yaml="a.yaml")
    actual = Path("composition_energies.yaml").read_text()
    expected = """H1:
  energy: 1.0
  source: mp-1
"""
    assert actual == expected
