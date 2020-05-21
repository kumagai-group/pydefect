# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from pydefect.corrections.manual_correction import NoCorrection, \
    ManualCorrection


def test_no_correction():
    assert NoCorrection().correction_energy == 0.0


def test_manual_correction():
    assert ManualCorrection(7.9).correction_energy == 7.9


"""
TODO
-

DONE
"""