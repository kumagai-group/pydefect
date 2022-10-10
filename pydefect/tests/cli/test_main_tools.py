# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import logging
from pathlib import Path

import pytest

from pydefect.cli.main_tools import sanitize_matrix, str_int_to_int, parse_dirs


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


def print_a_to_file_x(path: Path):
    filename = path / "x"
    filename.write_text("a")


def test_parse_dirs(tmpdir):
    pathlib_tmpdir = Path(tmpdir)
    parse_dirs([Path(tmpdir)], _inner_function=print_a_to_file_x)

    actual = pathlib_tmpdir / "x"
    assert actual.read_text() == "a"


def test_parse_dirs_file(tmpdir):
    tmpfile = Path(tmpdir) / "file"
    tmpfile.touch()
    parse_dirs([tmpfile], _inner_function=print_a_to_file_x)


def test_parse_dirs_fail(tmpdir, capsys):
    parse_dirs([Path(tmpdir) / "y"], _inner_function=print_a_to_file_x)
    captured = capsys.readouterr()
    assert captured.out == "No such file or directory\n"


def test_parse_dirs_skip_when_output_file_exist(tmpdir):
    LOGGER = logging.getLogger(__name__)
    LOGGER.info('Testing now.')
    (Path(tmpdir) / "x").touch()
    parse_dirs([Path(tmpdir)], _inner_function=print_a_to_file_x,
               output_filename="x")
