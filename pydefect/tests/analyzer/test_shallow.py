# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.analyzer.shallow import Shallow
from pydefect.tests.helpers.assertion import assert_msonable, \
    assert_json_roundtrip


@pytest.fixture
def shallow():
    return Shallow(is_shallow=True, method="manual")


def test_shallow_msonable(shallow):
    assert_msonable(shallow)


def test_shallow_to_json_file(shallow, tmpdir):
    assert_json_roundtrip(shallow, tmpdir)


"""
TODO
-

DONE
"""