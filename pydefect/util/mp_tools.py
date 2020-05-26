# -*- coding: utf-8 -*-

from pathlib import Path
from shutil import copyfile
from typing import List

import yaml
from pymatgen import Element, MPRester, Composition

from pydefect.chem_pot_diag.gas import MOLECULE_DATA

elements = [e.name for e in Element]
mol_dir = Path(__file__).parent / ".." / "chem_pot_diag" / "molecules"


class MpQuery:
    def __init__(self,
                 element_list: List[str],
                 properties: List[str] = None,
                 e_above_hull: float = 1e-5):

        self.element_list = element_list
        self.properties = (properties or
                           ["task_id", "full_formula", "final_energy",
                            "structure", "spacegroup", "band_gap",
                            "total_magnetization", "magnetic_type"])
        self.e_above_hull = e_above_hull

        excluded = list(set(elements) - set(element_list))
        # API key is parsed via .pmgrc.yaml

        with MPRester() as m:
            self.materials = \
                m.query(criteria={"elements": {"$in": element_list,
                                               "$nin": excluded},
                                  "e_above_hull": {"$lte": e_above_hull}},
                        properties=self.properties)


def make_poscars(materials_query: List[dict], path: Path = Path.cwd()) -> None:

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
            m.pop("structure").to(filename=m_path / "POSCAR")

            d = {"total_magnetization": m["total_magnetization"],
                 "band_gap": m["band_gap"],
                 "data_source": m["task_id"]}

            (m_path / "prior_info.yaml").write_text(yaml.dump(d))


