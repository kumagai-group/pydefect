# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from monty.json import MSONable
from monty.serialization import loadfn
from pymatgen import Element


@dataclass
class Energy(MSONable):
    name: str
    charge: int
    rel_energy: float
    atom_io: Dict[Element, int]
    correction_energy: Optional[Dict[str, float]] = None
    is_shallow: Optional[bool] = None

    def to_yaml(self, filename: str = "energy.yaml") -> None:
        lines = [f"name: {self.name}",
                 f"charge: {self.charge}",
                 f"rel_energy: {self.rel_energy}",
                 f"atom_io:"]
        for k, v in self.atom_io.items():
            lines.append(f"  {k}: {v}")
        if self.correction_energy is not None:
            lines.append(f"correction_energy:")
            for k, v in self.correction_energy.items():
                lines.append(f"  {k}: {v}")
        if self.is_shallow is not None:
            lines.append(f"is_shallow: {self.is_shallow}")
        Path(filename).write_text("\n".join(lines))

    @classmethod
    def from_yaml(cls, filename: str = "energy.yaml"):
        d = loadfn(filename)
        if d["atom_io"] is None:
            d["atom_io"] = {}
        else:
            d["atom_io"] = {Element(k): v for k, v in d["atom_io"].items()}
        return cls(**d)

    @property
    def total_correction(self):
        if self.correction_energy is not None:
            return sum([v for v in self.correction_energy.values()])
        return None

    def defect_formation_energy(self,
                                abs_chem_pot: Dict[Element, float]) -> float:
        if self.total_correction is None:
            raise ValueError("correction argument must be set to evaluate the"
                             "defect formation energy.")

        return (self.rel_energy + self.total_correction -
                reservoir_energy(self.atom_io, abs_chem_pot))


def reservoir_energy(atom_io: Dict[Element, int],
                     abs_chem_pot: Dict[Element, float]) -> float:
    return sum([-diff * abs_chem_pot[elem] for elem, diff in atom_io.items()])