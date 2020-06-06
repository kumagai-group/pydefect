# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pathlib import Path
from shutil import copyfile
from typing import List

import yaml
from pydefect.cli.vasp.molecules.molecules import MOLECULE_DATA
from pymatgen import Composition

mol_dir = Path(__file__).parent / "molecules"


def make_poscars_from_query(materials_query: List[dict], path: Path) -> None:

    for m in materials_query:
        reduced_formula = Composition(m["full_formula"]).reduced_formula

        if reduced_formula in MOLECULE_DATA.keys():
            dirname = path / f"mol_{reduced_formula}"
            if not dirname.exists():
                dirname.mkdir()
                copyfile(mol_dir / reduced_formula / "POSCAR",
                         dirname / "POSCAR")
                copyfile(mol_dir / reduced_formula / "prior_info.yaml",
                         dirname / "prior_info.yaml")

        else:
            m_path = path / f"{reduced_formula}_{m['task_id']}"
            m_path.mkdir()
            m["structure"].to(filename=m_path / "POSCAR")

            d = {"total_magnetization": m["total_magnetization"],
                 "band_gap": m["band_gap"],
                 "data_source": m["task_id"]}

            prior_info = m_path / "prior_info.yaml"
            prior_info.write_text(yaml.dump(d), None)
