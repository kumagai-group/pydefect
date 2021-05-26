# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from copy import copy

import pytest
from pydefect.analyzer.calc_summary import SingleCalcSummary, CalcSummary
from pydefect.analyzer.defect_structure_info import DefectType, SymmRelation
from vise.tests.helpers.assertion import assert_msonable


@pytest.fixture
def single_summary():
    return SingleCalcSummary(charge=0,
                             atom_io={"O": -1},
                             electronic_conv=True,
                             ionic_conv=True,
                             is_energy_strange=False,
                             same_config_from_init=True,
                             defect_type=str(DefectType.vacancy),
                             symm_relation=str(SymmRelation.same),
                             donor_phs=False,
                             acceptor_phs=False,
                             unoccupied_deep_state=True,
                             occupied_deep_state=True)


def test_single_calc_summary_msonable(single_summary):
    assert_msonable(single_summary)


def test_single_calc_summary_properties(single_summary):
    other = copy(single_summary)
    assert single_summary.same_atom_charge_io(other) is True
    other.charge = 1
    assert single_summary.same_atom_charge_io(other) is False

    assert single_summary.is_converged is True
    assert single_summary.is_proper_result is True
    assert single_summary.is_unusual is False


@pytest.fixture
def calc_summary(single_summary):
    second = SingleCalcSummary(charge=1,
                               atom_io={"Mg": -1},
                               electronic_conv=True,
                               ionic_conv=False,
                               is_energy_strange=False,
                               same_config_from_init=False,
                               defect_type=str(DefectType.vacancy_split),
                               symm_relation=str(SymmRelation.supergroup),
                               donor_phs=True,
                               acceptor_phs=False,
                               unoccupied_deep_state=True,
                               occupied_deep_state=True)

    return CalcSummary({"Va_O1_0": single_summary, "Va_Mg1_1": second})


def test_calc_summary_msonable(calc_summary):
    assert_msonable(calc_summary)


def test_calc_summary_str(calc_summary):
    print(calc_summary)
