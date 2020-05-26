# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass

from monty.json import MSONable
from monty.serialization import loadfn

from pydefect.util.mix_in import ToJsonFileMixIn


@dataclass
class TestTest(MSONable, ToJsonFileMixIn):
    a: str = "aaa"


def test_to_json_file_mix_in(tmpdir):
    actual = TestTest()
    tmpdir.chdir()
    actual.to_json_file()
    actual = loadfn("test_test.json")
    assert actual == TestTest()

"""
TODO
-

DONE
"""