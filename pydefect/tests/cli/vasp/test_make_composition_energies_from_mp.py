# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import pytest
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergy, \
    CompositionEnergies
from pydefect.cli.vasp.make_composition_energies_from_mp import \
    make_composition_energies_from_mp, remove_higher_energy_comp
from pymatgen.core import Composition, Element
from pymatgen.entries.computed_entries import ComputedEntry


@pytest.fixture
def mp_query_mock(mocker):
    mock = mocker.patch(
        "pydefect.cli.vasp.make_composition_energies_from_mp.MpEntries")
    mock.return_value.materials = \
        [ComputedEntry(Composition("O8"), -39.58364375, entry_id="mp-1"),
         ComputedEntry(Composition("Mg3"), -4.79068775, entry_id="mp-2"),
         ComputedEntry(Composition("Mg1O1"), -11.96742144, entry_id="mp-3")]
    return mock


@pytest.fixture
def comp_energies():
    return CompositionEnergies({
        Composition('O8'): CompositionEnergy(-39.58364375, "mp-1"),
        Composition('Mg3'): CompositionEnergy(-4.79068775, "mp-2"),
        Composition('Mg1O1'): CompositionEnergy(-11.96742144, "mp-3")})


# mp
# O:  -1.54794890
# Mg: -0.01731593
diff = {"Mg": -0.00912097 - -0.01731593, "O": -1.61154565 - -1.54794890}


@pytest.fixture
def comp_energies_corr():
    return CompositionEnergies({
        Composition('O8'): CompositionEnergy(-39.58364375 + diff["O"] * 8, "mp-1"),
        Composition('Mg3'): CompositionEnergy(-4.79068775 + diff["Mg"] * 3, "mp-2"),
        Composition('Mg1O1'): CompositionEnergy(-11.96742144 + diff["Mg"] + diff["O"], "mp-3")})


def test_make_chem_pot_diag_from_mp(mp_query_mock, comp_energies):
    actual = make_composition_energies_from_mp(elements=["Mg", "O"])
    mp_query_mock.assert_called_once_with(["Mg", "O"])
    assert actual == comp_energies


def test_make_chem_pot_diag_from_mp_yaml(mp_query_mock, tmpdir,
                                         comp_energies_corr):
    tmpdir.chdir()
    tmpdir.join("tmp.yaml").write("""O:  -1.61154565
Mg: -0.00912097""")
    actual = make_composition_energies_from_mp(elements=["Mg", "O"],
                                               atom_energy_yaml="tmp.yaml")
    assert actual == comp_energies_corr


def test_remove_higher_energy_comp():
    comp_es = CompositionEnergies(
        {Composition('Mg1'): CompositionEnergy(-1.9, "mp-2"),
         Composition('Mg2'): CompositionEnergy(-4.0, "mp-3"),
         Composition('O2'): CompositionEnergy(-4.0, "mp-4")})
    actual = remove_higher_energy_comp(comp_es)
    expected = CompositionEnergies(
        {Composition('Mg2'): CompositionEnergy(-4.0, "mp-3"),
         Composition('O2'): CompositionEnergy(-4.0, "mp-4"),})

    assert actual == expected

