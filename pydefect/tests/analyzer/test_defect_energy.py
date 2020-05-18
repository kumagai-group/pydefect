# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.analyzer.defect_energy import DefectEnergy


@pytest.fixture
def defect_energy():
    return DefectEnergy(name="Va_O1",
                        num_atom_change={"O": -1},
                        charges=[0, 1, 2],
                        energies=[6, 3, -4],
                        corrections=[0, 0, 0],
                        shallow=[False, False, False],
                        )


def test(defect_energy):
    actual = defect_energy.cross_points(1, 6)
    expected = ([[1.0, -2.0], [5.0, 6.0], [6.0, 6.0]], [[5.0, 6.0]])
    assert actual == expected


"""
TODO
- Evaluate the crossing points at given Fermi level range.

- Draw a single defect

DONE
"""