# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo, \
    NoElementEnergyError, CompositionEnergy, replace_comp_energy
# When same keys are inserted, only the latter one is accepted.
from pydefect.tests.helpers.assertion import assert_msonable
from pymatgen import Composition, Element

energies = [CompositionEnergy(Composition("H"), 0.0, "a"),
            CompositionEnergy(Composition("H4O2"), -4.0, "c"),
            CompositionEnergy(Composition("O"), 1.0, "b"),
            CompositionEnergy(Composition("Cl"), 12.0, "f"),
            CompositionEnergy(Composition("O2Cl"), 3.0, "e"),
            CompositionEnergy(Composition("O2Cl2"), 6.0, "d"),
            ]


@pytest.fixture
def cpd():
    return ChemPotDiag(energies, target=Composition("H2O"),
                       vertex_elements=[Element.H, Element.O])


def test_print(cpd):
    expected = """+----+--------+--------+---------+----------------------+
|    |   mu_H |   mu_O |   mu_Cl | Cl competing phase   |
|----+--------+--------+---------+----------------------|
| A  |    0   |     -3 |      -7 | ClO                  |
| B  |   -1.5 |      0 |     -11 | ClO2                 |
+----+--------+--------+---------+----------------------+"""
    assert cpd.__repr__() == expected


def test_composition_roundtrip(cpd):
    comp_e = CompositionEnergy(Composition("H"), 0.0, "a")
    rounded_comp_e = comp_e.__class__.from_dict(comp_e.as_dict())
    assert comp_e == rounded_comp_e


def test_msonable(cpd):
    assert_msonable(cpd)


def test_chem_pot_diag_yaml(cpd, tmpdir):
    tmpdir.chdir()
    cpd.to_yaml()
    expected = """Cl1:
  energy: 12.0
  source: f
H1:
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
vertex_elements:
- H
- O
"""
    assert Path("cpd.yaml").read_text() == expected

    actual = ChemPotDiag.from_yaml("cpd.yaml")
    assert actual.target == cpd.target


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
    assert cpd.vertex_elements_abs_energies_per_atom == {Composition("H2"): 0.0,
                                                         Composition("O2"): 1.0,
                                                         Composition("H2O"): -2/3}


def test_host_ele_abs_energies_per_atom_2():
    energies_with_cl = deepcopy(energies)
    energies_with_cl.append(CompositionEnergy(Composition("Cl"), 10.0, "z"))
    cpd = ChemPotDiag(energies_with_cl, target=Composition("H2O"),
                      vertex_elements=[Element.H, Element.O, Element.Cl])
    assert cpd.vertex_elements_abs_energies_per_atom == {Composition("H2"): 0.0,
                                                         Composition("O2"): 1.0,
                                                         Composition("H2O"): -2/3,
                                                         Composition("ClO"): 1.5,
                                                         Composition("ClO2"): 1.0,
                                                         Composition("Cl2"): 10.0}


def test_cpd_plot_info_lacking_element_data():
    new_energies = deepcopy(energies)
    new_energies.append(CompositionEnergy(Composition("MgO"), -3.0, "f"))
    with pytest.raises(NoElementEnergyError):
        ChemPotDiag(new_energies, target={"Mg": 1, "O": 1}).offset_to_abs


def test_chem_pot_diag_min_energy(cpd):
    assert cpd.lowest_relative_energy == -3


def test_impurity_abs_energy(cpd):
    expected = CompositionEnergy(Composition("Cl2O2"), 6.0, 'd'), 5.0
    assert cpd.impurity_abs_energy(Element.Cl, "A") == expected


def test_impurity_rel_energy(cpd):
    expected = CompositionEnergy(Composition("Cl2O2"), 6.0, 'd'), -7.0
    assert cpd.impurity_rel_energy(Element.Cl, "A") == expected


def test_vertex_list(cpd):
    index = ['A', 'B']
    columns = ['mu_H', 'mu_O', 'mu_Cl', 'Cl competing phase']
    data = [[0.0, -3.0, -7.0, 'ClO'], [-1.5, 0.0, -11.0, 'ClO2']]
    expected = pd.DataFrame(data, index=index, columns=columns)
    pd.testing.assert_frame_equal(cpd.target_vertex_list_dataframe, expected)


def test_target_formation_energy(cpd):
    assert cpd.target_formation_energy == ((-4.0) - (0.0 * 4 + 1.0 * 2)) / 6


@pytest.fixture
def cpd2():
    energies_here = [CompositionEnergy(Composition("H"), 0.0, "a"),
                     CompositionEnergy(Composition("H4O2"), -4.0, "c"),
                     CompositionEnergy(Composition("O"), 1.0, "b"),
                     CompositionEnergy(Composition("Cl2"), 3.0, "e"),
                     ]
    return ChemPotDiag(energies_here, target=Composition("H2O"),
                       vertex_elements=[Element.H, Element.O])


def test_impurity_abs_energy2(cpd2):
    expected = CompositionEnergy(Composition("Cl2"), 3.0, 'e'), 1.5
    assert cpd2.impurity_abs_energy(Element.Cl, "A") == expected


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


def test_replace_comp_energy():
    cpd = ChemPotDiag([CompositionEnergy(Composition("H"), 0.0, "a"),
                       CompositionEnergy(Composition("O"), 1.0, "b"),
                       CompositionEnergy(Composition("F"), 1.0, "c")],
                      target={"H": 1})
    replace_comp_energy(cpd,
                        [CompositionEnergy(Composition("H"), -1.0, "x"),
                         CompositionEnergy(Composition("O"), 0.0, "y")])

    expected = ChemPotDiag([CompositionEnergy(Composition("H"), -1.0, "x"),
                            CompositionEnergy(Composition("O"), 0.0, "y"),
                            CompositionEnergy(Composition("F"), 1.0, "c")],
                           target={"H": 1})
    assert cpd == expected


def test(tmpdir):
    tmpdir.chdir()
    tmpdir.join("cpd.yaml").write("""F8:
  energy: -6.69753983
  source: local
Mg2F4:
  energy: -31.54326603
  source: local
Mg3:
  energy: -4.50835875
  source: local
Na1F1:
  energy: -8.67062716
  source: local
Na20:
  energy: -25.85889772
  source: local
Na4Mg4F12:
  energy: -98.24594267
  source: local
target: MgF2
vertex_elements:
- Na
- Mg
- F
""")
    cpd = ChemPotDiag.from_yaml("cpd.yaml")
    expected = [[0, 0, -6.5404898],
                [-1.83646145e-01, 0, -6.35684365e+00],
                [-3.62484384e-01, 0, -6.29723090e+00],
                [-6.5404898, -12.7136873, 0],
                [-6.65971529, -12.59446181, 0]]
    actual = cpd.vertex_coords
    np.testing.assert_array_almost_equal(expected, actual)

"""
TODO
* Implement impurity_abs_energy
* implement print
"""