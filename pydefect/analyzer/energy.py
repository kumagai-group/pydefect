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
    formation_energy: float
    atom_io: Dict[Element, int]
    correction_energy: Dict[str, float]
    is_shallow: Optional[bool] = None

    @property
    def total_correction(self) -> float:
        return sum([v for v in self.correction_energy.values()])

    # def formation_energy_wo_corr(self,
    #                              abs_chem_pot: Dict[Element, float]) -> float:
    #     return self.rel_energy - reservoir_energy(self.atom_io, abs_chem_pot)


@dataclass
class EnergySummary(MSONable):
    name: str
    charge: int
    energy: Energy

    def to_yaml(self, filename: str = "energy_summary.yaml") -> None:
        lines = [f"name: {self.name}",
                 f"charge: {self.charge}",
                 f"formation_energy: {self.energy.formation_energy}",
                 f"atom_io:"]
        for k, v in self.energy.atom_io.items():
            lines.append(f"  {k}: {v}")
        lines.append(f"correction_energy:")
        for k, v in self.energy.correction_energy.items():
            lines.append(f"  {k}: {v}")
        is_shallow = "" if self.energy.is_shallow is None else self.energy.is_shallow
        lines.append(f"is_shallow: {is_shallow}")
        Path(filename).write_text("\n".join(lines))

    @classmethod
    def from_yaml(cls, filename: str = "energy.yaml") -> "EnergySummary":
        d = loadfn(filename)
        d["atom_io"] = {Element(k): v for k, v in d["atom_io"].items()} \
            if d["atom_io"] else {}
        return cls(name=d.pop("name"), charge=d.pop("charge"), energy=Energy(**d))

def reservoir_energy(atom_io: Dict[Element, int],
                     abs_chem_pot: Dict[Element, float]) -> float:
    return sum([-diff * abs_chem_pot[elem] for elem, diff in atom_io.items()])