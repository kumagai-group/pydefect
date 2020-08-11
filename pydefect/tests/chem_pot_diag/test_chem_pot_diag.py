# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo, \
    NoElementEnergyError, CompositionEnergy
# When same keys are inserted, only the latter one is accepted.
from pymatgen import Composition, Element

energies = {CompositionEnergy(Composition("H"), 0.0, "a"),
            CompositionEnergy(Composition("O"), 1.0, "b"),
            CompositionEnergy(Composition("H4O2"), -4.0, "c"),
            CompositionEnergy(Composition("O2Cl2"), 6.0, "d"),
            CompositionEnergy(Composition("O2Cl"), 3.0, "e")}


@pytest.fixture
def cpd():
    return ChemPotDiag(energies, target=Composition("H2O"))


def test_chem_pot_diag_yaml(cpd, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    cpd.to_yaml()
    expected = """H1:
  energy: 0.0
  source: a
H4O2:
  energy: -4.0
  source: c
O1:
  energy: 1.0
  source: b
O2Cl1:
  energy: 3.0
  source: e
O2Cl2:
  energy: 6.0
  source: d
target: H2O1
"""
    assert Path("cpd.yaml").read_text() == expected

    actual = ChemPotDiag.from_yaml("cpd.yaml")
    assert actual == cpd


def test_host_comp_abs_energies(cpd):
    assert cpd.rel_energies == {Composition("H2O1"): -1.0,
                                Composition("O2"): 0.0,
                                Composition("H2"): 0.0}


def test_impurity_elements(cpd):
    assert cpd.impurity_elements == [Element.Cl]


def test_chem_pot_diag(cpd):
    assert cpd.vertex_elements == [Element.H, Element.O]
    assert cpd.offset_to_abs == [0.0, 1.0]
    assert cpd.rel_energies[Composition("H2O")] == (-2.0 - 1.0) / 3
    np.testing.assert_array_almost_equal(cpd.vertex_coords,
                                         [[0.0, -3.0], [-1.5, 0.0]])
    np.testing.assert_array_almost_equal(cpd.target_vertices["A"], [0.0, -3.0])


def test_abs_chem_pots(cpd):
    assert cpd.abs_chem_pot_dict("A") == {Element.H: 0.0, Element.O: -3.0 + 1.0, Element.Cl: 5.0}


def test_host_ele_abs_energies_per_atom(cpd):
    assert cpd.host_ele_abs_energies_per_atom == {Composition("H2"): 0.0,
                                                  Composition("O2"): 1.0,
                                                  Composition("H2O"): -2/3}


def test_host_ele_abs_energies_per_atom_2():
    energies_with_cl = deepcopy(energies)
    energies_with_cl.add(CompositionEnergy(Composition("Cl"), 10.0, "z"))
    cpd = ChemPotDiag(energies_with_cl, target=Composition("H2O"),
                      host_elements=[Element.H, Element.O, Element.Cl])
    assert cpd.host_ele_abs_energies_per_atom == {Composition("H2"): 0.0,
                                                  Composition("O2"): 1.0,
                                                  Composition("H2O"): -2/3,
                                                  Composition("ClO"): 1.5,
                                                  Composition("ClO2"): 1.0,
                                                  Composition("Cl2"): 10.0}


def test_cpd_plot_info_lacking_element_data():
    new_energies = deepcopy(energies)
    new_energies.add(CompositionEnergy(Composition("MgO"), -3.0, "f"))
    with pytest.raises(NoElementEnergyError):
        ChemPotDiag(new_energies, target={"Mg": 1, "O": 1}).offset_to_abs


def test_chem_pot_diag_min_energy(cpd):
    assert cpd.lowest_relative_energy == -3


def test_impurity_abs_energy(cpd):
    expected = CompositionEnergy(Composition("Cl2O2"), 6.0, 'd'), 5.0
    assert cpd.impurity_abs_energy(Element.Cl, "A") == expected


@pytest.fixture()
def cpd_plot_info():
    cpd = ChemPotDiag(energies, target=Composition("H2O"))
    return CpdPlotInfo(cpd, min_range=-10)


@pytest.fixture
def cpd_plot_info_wo_min_range(cpd):
    return CpdPlotInfo(cpd)


def test_cpd_plot_info(cpd_plot_info):
    assert cpd_plot_info.comp_vertices == {Composition('H2'): [[0.0, -10.0], [0.0, -3.0]],
                                           Composition('H2O'): [[0.0, -3.0], [-1.5, 0.0]],
                                           Composition('O2'): [[-10.0, 0.0], [-1.5, 0.0]]}
    assert cpd_plot_info.comp_centers == {Composition('H2'): [0.0, -6.5],
                                          Composition('H2O'): [-0.75, -1.5],
                                          Composition('O2'): [-5.75, 0.0]}
    assert cpd_plot_info.atomic_fractions(Composition('H2O')) == [2 / 3, 1 / 3]


def test_cpd_plot_info_with_defaults(cpd_plot_info_wo_min_range):
    assert cpd_plot_info_wo_min_range.comp_vertices == \
           {Composition('H2'): [[0.0, -3.0], [0.0, -3.3]],
            Composition('H2O'): [[-1.5, 0.0], [0.0, -3.0]],
            Composition('O2'): [[-3.3, 0.0], [-1.5, 0.0]]}


"""
TODO
1. Implement impurity_abs_energy
4. implement print


DONE
2. change energies types from Dict[str, float] to CompoundEnergy
3. to yaml
"""