# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import copy

import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, \
    CompositionEnergy, CompositionEnergies, StandardEnergies, RelativeEnergies, \
    ChemPotDiagMaker, TargetVertices, TargetVertex, \
    ImpurityInfo, calc_chem_pot
# When same keys are inserted, only the latter one is accepted.
from pymatgen.core import Composition
from vise.tests.helpers.assertion import assert_yaml_roundtrip


@pytest.fixture
def comp_es():
    return CompositionEnergies(
        {Composition("H"): CompositionEnergy(0.0, "a"),
         Composition("H4O2"): CompositionEnergy(-4.0, "c"),
         Composition("O"): CompositionEnergy(1.0, "b"),
         Composition("Cl"): CompositionEnergy(12.0, "f"),
         Composition("O2Cl"): CompositionEnergy(3.0, "e"),
         Composition("O2Cl2"): CompositionEnergy(6.0, "d")})


def test_comp_es_yaml_roundtrip(comp_es, tmpdir):
    expected_text = """Cl1:
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
"""
    assert_yaml_roundtrip(comp_es, tmpdir, expected_text, pass_comp_dict=True, pass_comp_items=True)


def test_comp_es_elements(comp_es):
    assert comp_es.elements == ["Cl", "H", "O"]


def test_std_rel_energies(comp_es):
    actual = comp_es.std_rel_energies
    expected_ref = StandardEnergies({"H": 0.0, "O": 1.0, "Cl": 12.0})
    expected_rel = RelativeEnergies({'ClO': -5.0, 'ClO2': -3.666666666666667, 'H2O': -1.0})
    assert actual[0] == expected_ref
    assert actual[1] == expected_rel


@pytest.fixture
def std_energies():
    return StandardEnergies({"Mg": 1.0, "O": 2.0})


def test_std_energies(std_energies, tmpdir):
    expected_text = """Mg: 1.0
O: 2.0
"""
    assert_yaml_roundtrip(std_energies, tmpdir, expected_text, pass_comp_dict=True, pass_comp_items=True)


def test_calc_chem_pot():
    actual = calc_chem_pot("MgO3", energy_per_atom=-10.0, target_elem="Mg",
                           other_elem_chem_pot={"O": -4.0})
    assert actual == (-10.0 - (-4.0 * 3 / 4)) * 4


@pytest.fixture
def relative_energies():
    return RelativeEnergies({"MgO2": -3.0, "MgAlO2": -1.0})


def test_relative_energies_yaml(relative_energies, tmpdir):
    expected_text = """MgAlO2: -1.0
MgO2: -3.0
"""
    assert_yaml_roundtrip(relative_energies, tmpdir, expected_text, pass_comp_dict=True, pass_comp_items=True)


def test_relative_energies_element_set(relative_energies):
    assert relative_energies.all_element_set == {"Mg", "Al", "O"}
    # assert relative_energies.impurity_element_set(["Mg", "O"]) == {"O"}


def test_relative_energies_related_comp_energies(relative_energies):
    actual = relative_energies.related_comp_energies(elements=["Mg", "O"])
    expected = {"MgO2": -3.0}
    assert actual == expected


def test_relative_energies_comp_energies_w_element(relative_energies):
    actual = relative_energies.comp_energies_w_element(element="Al")
    expected = {"MgAlO2": -1.0}
    assert actual == expected


@pytest.fixture
def cpd_maker():
    return ChemPotDiagMaker(RelativeEnergies({"MgO2": -3.0, "MgAlO2": -2.5}),
                            elements=["Mg", "O"])


def test_cpd_maker_chem_pot_diag(cpd_maker):
    actual = cpd_maker.chem_pot_diag
    min_val = -9.0 * 1.1
    expected = ChemPotDiag(vertex_elements=["Mg", "O"],
                           polygons={'Mg': [[0.0, min_val], [0.0, -4.5]],
                                     'MgO2': [[0.0, -4.5], [-9.0, 0.0]],
                                     'O': [[min_val, 0.0], [-9.0, 0.0]]})
    assert actual == expected


@pytest.fixture
def cpd():
    vertex_a = TargetVertex(coords=[0.0, -4.5],
                            competing_phases=["Mg"],
                            impurity_info={"Al": ImpurityInfo(-1.0, "MgAlO2")})
    vertex_b = TargetVertex(coords=[-9.0, 0.0],
                            competing_phases=["O"],
                            impurity_info={"Al": ImpurityInfo(-1.0, "MgAlO2")})
    min_val = -9.0 * 1.1
    return ChemPotDiag(
        vertex_elements=["Mg", "O"],
        polygons={'Mg': [[0.0, min_val], [0.0, -4.5]],
                  'MgO2': [[0.0, -4.5], [-9.0, 0.0]],
                  'O': [[min_val, 0.0], [-9.0, 0.0]]},
        target_vertices=TargetVertices("MgO2", {"A": vertex_a, "B": vertex_b}))


def test_cpd_maker_chem_pot_diag_with_target(cpd_maker, cpd):
    cpd_maker_2 = copy(cpd_maker)
    cpd_maker_2.target = "MgO2"
    actual = cpd_maker_2.chem_pot_diag
    assert actual == cpd


def test_cpd(cpd):
    assert cpd.dim == 2
    assert cpd.comp_centers == {"Mg": [0.0, (-9.9 - 4.5) / 2],
                                "MgO2": [-4.5, -2.25],
                                "O": [(-9.9 - 9.0) / 2, 0.0]}
    assert cpd.min_value == -9.9

# def test_print(cpd_h2o):
#     expected = """+----+--------+--------+---------+----------------------+
# |    |   mu_H |   mu_O |   mu_Cl | Cl competing phase   |
# |----+--------+--------+---------+----------------------|
# | A  |    0   |     -3 |      -7 | ClO                  |
# | B  |   -1.5 |      0 |     -11 | ClO2                 |
# +----+--------+--------+---------+----------------------+"""
#     assert cpd_h2o.__repr__() == expected
#
#
#
#
#
# def test_composition_roundtrip(cpd_h2o):
#     comp_e = CompositionEnergy(Composition("H"), 0.0, "a")
#     rounded_comp_e = comp_e.__class__.from_dict(comp_e.as_dict())
#     assert comp_e == rounded_comp_e
#
#
# def test_msonable(cpd_h2o):
#     assert_msonable(cpd_h2o)
#
#
# def test_chem_pot_diag_yaml(cpd_h2o, tmpdir):
#     tmpdir.chdir()
#     cpd_h2o.to_yaml()
#     expected = """Cl1:
#   energy: 12.0
#   source: f
# H1:
#   energy: 0.0
#   source: a
# H4O2:
#   energy: -4.0
#   source: c
# O1:
#   energy: 1.0
#   source: b
# O2Cl1:
#   energy: 3.0
#   source: e
# O2Cl2:
#   energy: 6.0
#   source: d
# target: H2O1
# vertex_elements:
# - H
# - O
# """
#     assert Path("cpd.yaml").read_text() == expected
#
#     actual = ChemPotDiag.from_yaml("cpd.yaml")
#     assert actual.target == cpd_h2o.target
#
#
# def test_host_comp_abs_energies(cpd_h2o):
#     assert cpd_h2o.rel_energies == {Composition("H2O1"): -1.0,
#                                 Composition("O2"): 0.0,
#                                 Composition("H2"): 0.0}
#
#
# def test_impurity_elements(cpd_h2o):
#     assert cpd_h2o.impurity_elements == [Element.Cl]
#
#
# def test_chem_pot_diag(cpd_h2o):
#     assert cpd_h2o.vertex_elements == [Element.H, Element.O]
#     assert cpd_h2o.offset_to_abs == [0.0, 1.0]
#     assert cpd_h2o.rel_energies[Composition("H2O")] == (-2.0 - 1.0) / 3
#     np.testing.assert_array_almost_equal(cpd_h2o.vertex_coords,
#                                          [[0.0, -3.0], [-1.5, 0.0]])
#     np.testing.assert_array_almost_equal(cpd_h2o.target_vertices["A"], [0.0, -3.0])
#
#
# def test_abs_chem_pots(cpd_h2o):
#     assert cpd_h2o.abs_chem_pot_dict("A") == {Element.H: 0.0, Element.O: -3.0 + 1.0, Element.Cl: 5.0}
#     assert cpd_h2o.label_at_rich_condition(Element.H) == "A"
#     assert cpd_h2o.label_at_rich_condition(Element.O) == "B"
#     assert cpd_h2o.label_at_poor_condition(Element.O) == "A"
#
#
# def test_host_ele_abs_energies_per_atom(cpd_h2o):
#     assert cpd_h2o.vertex_elements_abs_energies_per_atom == {Composition("H2"): 0.0,
#                                                          Composition("O2"): 1.0,
#                                                          Composition("H2O"): -2/3}
#
#
# def test_host_ele_abs_energies_per_atom_2():
#     energies_with_cl = deepcopy(energies)
#     energies_with_cl.append(CompositionEnergy(Composition("Cl"), 10.0, "z"))
#     cpd = ChemPotDiag(energies_with_cl, target=Composition("H2O"),
#                       vertex_elements=[Element.H, Element.O, Element.Cl])
#     assert cpd.vertex_elements_abs_energies_per_atom == {Composition("H2"): 0.0,
#                                                          Composition("O2"): 1.0,
#                                                          Composition("H2O"): -2/3,
#                                                          Composition("ClO"): 1.5,
#                                                          Composition("ClO2"): 1.0,
#                                                          Composition("Cl2"): 10.0}
#
#
# def test_cpd_plot_info_lacking_element_data():
#     new_energies = deepcopy(energies)
#     new_energies.append(CompositionEnergy(Composition("MgO"), -3.0, "f"))
#     with pytest.raises(NoElementEnergyError):
#         ChemPotDiag(new_energies, target={"Mg": 1, "O": 1}).offset_to_abs
#
#
# def test_chem_pot_diag_min_energy(cpd_h2o):
#     assert cpd_h2o.lowest_relative_energy == -3
#
#
# def test_impurity_abs_energy(cpd_h2o):
#     expected = CompositionEnergy(Composition("Cl2O2"), 6.0, 'd'), 5.0
#     assert cpd_h2o.impurity_abs_energy(Element.Cl, "A") == expected
#
#
# def test_impurity_rel_energy(cpd_h2o):
#     expected = CompositionEnergy(Composition("Cl2O2"), 6.0, 'd'), -7.0
#     assert cpd_h2o.impurity_rel_energy(Element.Cl, "A") == expected
#
#
# def test_vertex_list(cpd_h2o):
#     index = ['A', 'B']
#     columns = ['mu_H', 'mu_O', 'mu_Cl', 'Cl competing phase']
#     data = [[0.0, -3.0, -7.0, 'ClO'], [-1.5, 0.0, -11.0, 'ClO2']]
#     expected = pd.DataFrame(data, index=index, columns=columns)
#     pd.testing.assert_frame_equal(cpd_h2o.target_vertex_list_dataframe, expected)
#
#
# def test_target_formation_energy(cpd_h2o):
#     assert cpd_h2o.target_formation_energy == ((-4.0) - (0.0 * 4 + 1.0 * 2)) / 6
#
#
# @pytest.fixture
# def cpd2():
#     energies_here = [CompositionEnergy(Composition("H"), 0.0, "a"),
#                      CompositionEnergy(Composition("H4O2"), -4.0, "c"),
#                      CompositionEnergy(Composition("O"), 1.0, "b"),
#                      CompositionEnergy(Composition("Cl2"), 3.0, "e"),
#                      ]
#     return ChemPotDiag(energies_here, target=Composition("H2O"),
#                        vertex_elements=[Element.H, Element.O])
#
#
# def test_impurity_abs_energy2(cpd2):
#     expected = CompositionEnergy(Composition("Cl2"), 3.0, 'e'), 1.5
#     assert cpd2.impurity_abs_energy(Element.Cl, "A") == expected
#
#
# @pytest.fixture
# def cpd_all_vertices():
#     energies_here = [CompositionEnergy(Composition("H"), 0.0, "a"),
#                      CompositionEnergy(Composition("H4O2"), -4.0, "c"),
#                      CompositionEnergy(Composition("O"), 1.0, "b"),
#                      CompositionEnergy(Composition("Cl2"), 3.0, "e"),
#                      ]
#     return ChemPotDiag(energies_here, target=Composition("H2O"),
#                        vertex_elements=[Element.H, Element.O, Element.Cl])
#
#
# def test_impurity_abs_energy3(cpd_all_vertices):
#     assert cpd_all_vertices.impurity_elements == []
#     with pytest.raises(ValueError):
#         cpd_all_vertices.impurity_abs_energy(Element.Cl, "A")
#     expected = {Element.H: -1.4551915228366852e-11,
#                 Element.O: -2.0,
#                 Element.Cl: 1.499999999985448}
#     assert cpd_all_vertices.abs_chem_pot_dict("A") == expected
#
#
# @pytest.fixture()
# def cpd_plot_info():
#     cpd = ChemPotDiag(energies, target=Composition("H2O"))
#     return CpdPlotInfo(cpd, min_range=-10)
#
#
# @pytest.fixture
# def cpd_plot_info_wo_min_range(cpd_h2o):
#     return CpdPlotInfo(cpd_h2o)
#
#
# def test_cpd_plot_info(cpd_plot_info):
#     assert cpd_plot_info.comp_vertices == {Composition('H2'): [[0.0, -10.0], [0.0, -3.0]],
#                                            Composition('H2O'): [[0.0, -3.0], [-1.5, 0.0]],
#                                            Composition('O2'): [[-10.0, 0.0], [-1.5, 0.0]]}
#     assert cpd_plot_info.comp_centers == {Composition('H2'): [0.0, -6.5],
#                                           Composition('H2O'): [-0.75, -1.5],
#                                           Composition('O2'): [-5.75, 0.0]}
#     assert cpd_plot_info.atomic_fractions(Composition('H2O')) == [2 / 3, 1 / 3]
#
#
# def test_cpd_plot_info_with_defaults(cpd_plot_info_wo_min_range):
#     assert cpd_plot_info_wo_min_range.comp_vertices == \
#            {Composition('H2'): [[0.0, -3.0], [0.0, -3.3]],
#             Composition('H2O'): [[-1.5, 0.0], [0.0, -3.0]],
#             Composition('O2'): [[-3.3, 0.0], [-1.5, 0.0]]}
#
#
# def test_replace_comp_energy():
#     cpd = ChemPotDiag([CompositionEnergy(Composition("H"), 0.0, "a"),
#                        CompositionEnergy(Composition("O"), 1.0, "b"),
#                        CompositionEnergy(Composition("F"), 1.0, "c")],
#                       target={"H": 1})
#     replace_comp_energy(cpd,
#                         [CompositionEnergy(Composition("H"), -1.0, "x"),
#                          CompositionEnergy(Composition("O"), 0.0, "y")])
#
#     expected = ChemPotDiag([CompositionEnergy(Composition("H"), -1.0, "x"),
#                             CompositionEnergy(Composition("O"), 0.0, "y"),
#                             CompositionEnergy(Composition("F"), 1.0, "c")],
#                            target={"H": 1})
#     assert cpd == expected
#
#
# def test(tmpdir):
#     tmpdir.chdir()
#     tmpdir.join("cpd.yaml").write("""F8:
#   energy: -6.69753983
#   source: local
# Mg2F4:
#   energy: -31.54326603
#   source: local
# Mg3:
#   energy: -4.50835875
#   source: local
# # Omit 1
# NaF:
#   energy: -8.67062716
#   source: local
# Na20:
#   energy: -25.85889772
#   source: local
# # Omit source
# Na4Mg4F12:
#   energy: -98.24594267
# target: MgF2
# vertex_elements:
# - Na
# - Mg
# - F
# """)
#     cpd = ChemPotDiag.from_yaml("cpd.yaml")
#     expected = [[0, 0, -6.5404898],
#                 [-1.83646145e-01, 0, -6.35684365e+00],
#                 [-3.62484384e-01, 0, -6.29723090e+00],
#                 [-6.5404898, -12.7136873, 0],
#                 [-6.65971529, -12.59446181, 0]]
#     actual = cpd.vertex_coords
#     np.testing.assert_array_almost_equal(expected, actual)
#
# """
# TODO
# * Implement impurity_abs_energy
# * implement print
# """