# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy

import numpy as np
import pytest
from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo, \
    NoElementEnergyError
# When same keys are inserted, only the latter one is accepted.
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip
from pymatgen import Composition, Element

energies = {"H": 0.0, "O": 1.0, "H4O2": -4.0}


@pytest.fixture
def cpd():
    return ChemPotDiag(energies, target=Composition("H2"))


def test_chem_pot_diag_assert_msonable(cpd):
    assert_msonable(cpd)


# 2020.5.30: Composition is not recovered with MSONable
def test_msonable_for_composition(cpd):
    a = ChemPotDiag.from_dict(cpd.as_dict())
    assert isinstance(a.target, Composition)


def test_chem_pot_diag_to_json_file(cpd, tmpdir):
    assert_json_roundtrip(cpd, tmpdir)


def test_chem_pot_diag(cpd):
    assert cpd.vertex_elements == [Element.H, Element.O]
    assert cpd.offset_to_abs == [0.0, 1.0]
    assert cpd.rel_energies[Composition("H2O")] == (-2.0 - 1.0) / 3
    np.testing.assert_array_almost_equal(cpd.vertex_coords,
                                         [[0.0, -3.0], [-1.5, 0.0]])
    np.testing.assert_array_almost_equal(cpd.target_vertices["A"], [0.0, -3.0])


def test_abs_chem_pots(cpd):
    assert cpd.abs_chem_pot_dict("A") == {Element.H: 0.0, Element.O: -3.0 + 1.0}


def test_cpd_plot_info_lacking_element_data():
    new_energies = deepcopy(energies)
    new_energies[Composition("MgO")] = -3.0
    with pytest.raises(NoElementEnergyError):
        ChemPotDiag(new_energies, target={"Mg": 1, "O": 1}).offset_to_abs


def test_chem_pot_diag_min_energy(cpd):
    assert cpd.lowest_relative_energy == -3


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
- Return the formation energy of a compound
- Return the absolute energy of a compound
- 

DONE
"""