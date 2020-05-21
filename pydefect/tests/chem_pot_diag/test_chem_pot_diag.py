# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from copy import deepcopy

import numpy as np
import pytest
from pymatgen import Composition, Element

from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo, \
    NoElementEnergyError

# When same keys are inserted, only the latter one is accepted.
energies = {Composition("H"): 0.0,
            Composition("O"): 0.5,
            Composition("O"): 1.0,
            Composition("H4O2"): -4.0}


@pytest.fixture
def cpd():
    return ChemPotDiag(energies, target=Composition("H2"))


def test_chem_pot_diag(cpd):
    assert cpd.vertex_elements == [Element.H, Element.O]
    assert cpd.offset_to_abs == [0.0, 1.0]
    assert cpd.rel_energies[Composition("H2O")] == (-2.0 - 1.0) / 3
    np.testing.assert_array_almost_equal(cpd.vertex_coords,
                                         [[0.0, -3.0], [-1.5, 0.0]])
    np.testing.assert_array_almost_equal(cpd.target_vertices["A"], [0.0, -3.0])


def test_cpd_plot_info_lacking_element_data():
    new_energies = deepcopy(energies)
    new_energies[Composition("MgO")] = -3.0
    with pytest.raises(NoElementEnergyError):
        ChemPotDiag(new_energies)


def test_chem_pot_diag_min_energy(cpd):
    assert cpd.min_rel_energies == -3.3


@pytest.fixture(scope="session")
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