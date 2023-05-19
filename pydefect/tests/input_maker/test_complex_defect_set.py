# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from pathlib import Path

import pytest
from pydefect.input_maker.complex_defect_set import ComplexDefectSet, \
    ComplexDefect


@pytest.fixture
def complex_defect():
    return ComplexDefect("comp", {1: None, 2: "N"}, charges=[1, 2])


@pytest.fixture
def complex_defect_set(complex_defect):
    return ComplexDefectSet(defects=[complex_defect])


def test_complex_defect_set(complex_defect_set):
    # test iterable
    for defect in complex_defect_set:
        defect


def test_from_yaml(tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    Path("tmp.yaml").write_text("""comp:
  replaced_sites:
    1: null
    2: N
  charges: 
    - 1
    - 2""")
    actual = ComplexDefectSet.from_yaml("tmp.yaml")
    assert actual.defects[0].name == "comp"
