# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from pathlib import Path

import pytest
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.concentration.degeneracy import Degeneracies, Degeneracy, \
    MakeDegeneracy
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.defect_structure_info import DefectStructureInfo

degeneracies = Degeneracies({"Va_O1": {0: Degeneracy(1, 1, "1", "1")}})


def test_degeneracies():
    assert degeneracies["Va_O1"] == {0: Degeneracy(1, 1, "1", "1")}


def test_degeneracies_yaml_round_trip(tmpdir):
    tmpdir.chdir()
    text = """Va_O1:
  0:
    final_site_sym: '1'
    initial_site_sym: '1'
    site: 1
    spin: 1
"""
    Path("degeneracies.yaml").write_text(text)
    actual = Degeneracies.from_yaml("degeneracies.yaml")
    assert isinstance(actual["Va_O1"][0], Degeneracy)


def test_make_degeneracy(mocker):
    make_deg = MakeDegeneracy("Pm-3m")
    assert make_deg.mag_to_spin_degeneracy(mag=0.09) == 1
    assert make_deg.mag_to_spin_degeneracy(mag=-1.0) == 3
    with pytest.raises(ValueError):
        make_deg.mag_to_spin_degeneracy(mag=0.89)

    m_energy_info = mocker.Mock(DefectEnergyInfo)
    m_energy_info.name, m_energy_info.charge = "Va_O1", 0

    m_calc_results = mocker.Mock(CalcResults)
    m_calc_results.magnetization = 1.0

    m_structure_info = mocker.Mock(DefectStructureInfo)
    m_structure_info.initial_site_sym = "4mmm"
    m_structure_info.final_site_sym = "2mm"

    assert make_deg.degeneracies == Degeneracies({})
    make_deg.add_degeneracy(m_energy_info, m_calc_results, m_structure_info)
    assert make_deg.degeneracies["Va_O1"][0] == Degeneracy(12, 3, "4mmm", "2mm")


