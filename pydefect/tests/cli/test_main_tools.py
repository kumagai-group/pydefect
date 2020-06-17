# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import pytest

from pydefect.cli.main_tools import sanitize_matrix, str_int_to_int


def test_sanitize_matrix_9_input_values():
    actual = sanitize_matrix(list(range(9)))
    expected = [[0, 1, 2], [3, 4, 5], [6, 7, 8]]
    assert actual == expected


def test_sanitize_matrix_3_input_values():
    actual = sanitize_matrix(list(range(1, 4)))
    expected = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
    assert actual == expected


def test_sanitize_matrix_1_input_value():
    actual = sanitize_matrix([4])
    expected = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
    assert actual == expected


def test_sanitize_matrix_raise_error_for_improper_matrix_element_length():
    with pytest.raises(ValueError):
         sanitize_matrix([1, 4])


def test_str_int_to_int():
    assert str_int_to_int("a") == "a"
    assert str_int_to_int("1") == 1
    assert str_int_to_int("1.1") == "1.1"


"""
TODO
- return matrix

DONE
"""