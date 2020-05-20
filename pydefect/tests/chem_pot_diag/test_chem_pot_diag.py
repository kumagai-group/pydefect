# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pymatgen import Composition

from pydefect.chem_pot_diag.chem_pot_diag import ChemPotDiag, CpdPlotInfo

# When same keys are inserted, only the latter one is accepted.
energies = {Composition("H"): 0.0,
            Composition("O"): 0.5,
            Composition("O"): 1.0,
            Composition("H2O"): -2.0}


@pytest.fixture
def cpd():
    return ChemPotDiag(energies)


def test_chem_pot_diag(cpd):
    assert cpd.vertex_elements == ["H", "O"]
    assert cpd.offset_to_abs == [0.0, 1.0]
    assert cpd.rel_energies[Composition("H2O")] == (-2.0 - 1.0) / 3
    np.testing.assert_array_almost_equal(cpd.vertex_coords,
                                         [[0.0, -3.0], [-1.5, 0.0]])


def test_chem_pot_diag_min_energy(cpd):
    assert cpd.min_rel_energies == -3.3


def test_cpd_plot_info(cpd):
    cpd_plot_info = CpdPlotInfo(cpd, target=Composition("H2O"), min_range=-10)
    assert cpd_plot_info.comp_vertices == {'H2': [[0.0, -10.0], [0.0, -3.0]],
                                           'H2O': [[0.0, -3.0], [-1.5, 0.0]],
                                           'O2': [[-10.0, 0.0], [-1.5, 0.0]]}
    assert cpd_plot_info.target_vertices == {"A": [0.0, -3.0], "B": [-1.5, 0.0]}


def test_cpd_plot_info_with_defaults(cpd):
    cpd_plot_info = CpdPlotInfo(cpd)
    assert cpd_plot_info.comp_vertices == {'H2': [[0.0, -3.0], [0.0, -3.3]],
                                           'H2O': [[-1.5, 0.0], [0.0, -3.0]],
                                           'O2': [[-3.3, 0.0], [-1.5, 0.0]]}

"""
TODO
- Return the formation energy of a compound
- Return the absolute energy of a compound
- 

DONE
"""