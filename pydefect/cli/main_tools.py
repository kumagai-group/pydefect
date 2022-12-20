# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import traceback
from pathlib import Path
from typing import List, Callable, Any

import numpy as np
from vise.util.logger import get_logger

logger = get_logger(__name__)


def sanitize_matrix(matrix: List[int]) -> List[List[int]]:

    if len(matrix) == 9:
        return [matrix[:3], matrix[3:6], matrix[6:]]
    elif len(matrix) == 3:
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[i]
        return result.tolist()
    elif len(matrix) == 1:
        result = np.eye(3, dtype=int)
        for i in range(3):
            result[i, i] = matrix[0]
        return result.tolist()
    else:
        raise ValueError(f"Matrix element length {len(matrix)} is improper.")


def str_int_to_int(x):
    try:
        return int(x)
    except ValueError:
        return x


def parse_dirs(dirs: List[Path],
               _inner_function: Callable[[Path], Any],
               verbose: bool = False,
               output_filename: str = None):
    failed_directories = []
    parsed_results = []
    for _dir in dirs:
        if _dir.is_file():
            logger.info(f"{_dir} is a file, so skipped.")
            continue
        if output_filename and (_dir / output_filename).exists():
            logger.info(f"In {_dir}, {output_filename} already exists.")
            continue

        logger.info(f"Parsing data in {_dir} ...")
        try:
            _return = _inner_function(_dir)
            if _return:
                parsed_results.append(_return)
        except Exception as e:
            if verbose:
                print(traceback.print_exc())
            else:
                print(e.args[1])
            logger.warning(f"Failing parsing {_dir} ...")
            failed_directories.append(str(_dir))
            continue
    if failed_directories:
        failed_dir_string = '\n'.join(failed_directories)
        logger.warning(f"Failed directories are:\n{failed_dir_string}")
        if verbose is False:
            logger.warning(f"To see details, try to run with --verbose")
    else:
        logger.info("Parsing all the directories succeeded.")

    return parsed_results or None