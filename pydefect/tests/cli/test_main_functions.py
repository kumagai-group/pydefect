# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    CompositionEnergy, RelativeEnergies, ChemPotDiag
from pydefect.cli.main_functions import make_std_and_rel_energies, \
    make_cpd_and_vertices, plot_chem_pot_diag
from pymatgen.core import Composition


def test_make_std_and_rel_energies(mocker, tmpdir):
    tmpdir.chdir()
    print(tmpdir)
    args = Namespace(composition_energies_yaml="a.yaml")
    mock = mocker.patch("pydefect.cli.main_functions.CompositionEnergies")
    mock.from_yaml.return_value = CompositionEnergies(
        {Composition("H"): CompositionEnergy(-1.0),
         Composition("He"): CompositionEnergy(-2.0),
         Composition("HHe"): CompositionEnergy(-5.0)})
    make_std_and_rel_energies(args)
    mock.from_yaml.assert_called_once_with("a.yaml")
    actual = Path("standard_energies.yaml").read_text()
    expected = """H: -1.0
He: -2.0
"""
    assert actual == expected

    actual = Path("relative_energies.yaml").read_text()
    expected = """HHe: -1.0
"""
    assert actual == expected


def test_make_cpd_and_vertices(mocker, tmpdir):
    tmpdir.chdir()
    print(tmpdir)
    mock = mocker.patch("pydefect.cli.main_functions.RelativeEnergies")
    mock.from_yaml.return_value = RelativeEnergies({"HHe": -1.0})

    args = Namespace(rel_energy_yaml="a.yaml", target="HHe", elements=None)
    make_cpd_and_vertices(args)

    actual = Path("target_vertices.yaml").read_text()
    expected = """target: HHe
A:
  chem_pot:
    H: 0.0
    He: -2.0
  competing_phases:
  - H
  impurity_phases: []
B:
  chem_pot:
    H: -2.0
    He: 0.0
  competing_phases:
  - He
  impurity_phases: []
"""
    assert actual == expected


def test_plot_chem_pot_diag(tmpdir):
    tmpdir.chdir()
    args = Namespace(chem_pot_diag=ChemPotDiag(vertex_elements=["Mg", "O"],
                                               polygons={"Mg": [[0.0, -3.0], [0.0, -1.0]],
                                                         "O": [[-3.0, 0.0], [-1.0, 0.0]],
                                                         "MgO": [[-1.0, 0.0], [0.0, -1.0]]}))
    plot_chem_pot_diag(args)
