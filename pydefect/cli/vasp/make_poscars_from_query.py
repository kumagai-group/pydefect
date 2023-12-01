# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path
from shutil import copyfile
from typing import List, Union

import yaml
from pydefect.cli.vasp.molecules.molecules import MOLECULE_DATA
from pymatgen.core import Composition
from vise.util.logger import get_logger

mol_dir = Path(__file__).parent / "molecules"

logger = get_logger(__name__)


def make_poscars_from_query(materials_query: List[Union[dict, "SummaryDoc"]],
                            path: Path) -> None:
    for query in materials_query:
        try:
            reduced_formula = Composition(query["full_formula"]).reduced_formula
        except TypeError:
            # For new MPRester.
            reduced_formula = query.formula_pretty

        if reduced_formula in MOLECULE_DATA.keys():
            _make_molecular_directory(path, reduced_formula)
        else:
            _make_solid_directory(path, reduced_formula, query)


def _make_solid_directory(path: Path, reduced_formula: str,
                          query: Union[dict, "SummaryDoc"]):
    try:
        task_id = query['task_id']
    except TypeError:
        # For new MPRester.
        task_id = query.material_id

    dirname = path / f"{reduced_formula}_{task_id}"
    logger.info(f"{dirname.name} is being created with a MP structure.")
    dirname.mkdir()

    try:
        query["structure"].to(filename=str(dirname / "POSCAR"))
        d = {"total_magnetization": query["total_magnetization"],
             "band_gap": query["band_gap"],
             "data_source": query["task_id"]}
    except TypeError:
        # For new MPRester.
        query.structure.to(filename=str(dirname / "POSCAR"))
        d = {"total_magnetization": query.total_magnetization,
             "band_gap": query.band_gap,
             "data_source": f"new MPRester {query.material_id}"}

    prior_info = dirname / "prior_info.yaml"
    prior_info.write_text(yaml.dump(d), None)


def _make_molecular_directory(path: Path, reduced_formula: str):
    dirname = path / f"mol_{reduced_formula}"
    logger.info(
        f"{dirname.name} is being created with molecular data in pydefect.")
    if not dirname.exists():
        dirname.mkdir()
        copyfile(mol_dir / reduced_formula / "POSCAR",
                 dirname / "POSCAR")
        copyfile(mol_dir / reduced_formula / "prior_info.yaml",
                 dirname / "prior_info.yaml")
