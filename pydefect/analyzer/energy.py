# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from monty.json import MSONable
from monty.serialization import loadfn


@dataclass
class Energy(MSONable):
    rel_energy: float
    atom_io: Dict[str, int]
    correction_energy: Optional[float] = None
    is_shallow: Optional[bool] = None

    def to_yaml(self, filename: str = "energy.yaml") -> None:
        lines = [f"rel_energy: {self.rel_energy}",
                 f"atom_io:"]
        for k, v in self.atom_io.items():
            lines.append(f"  {k}: {v}")
        if self.correction_energy is not None:
            lines.append(f"correction_energy: {self.correction_energy}")
        if self.is_shallow is not None:
            lines.append(f"is_shallow: {self.is_shallow}")
        Path(filename).write_text("\n".join(lines))

    @classmethod
    def from_yaml(cls, filename: str = "energy.yaml"):
        d = loadfn(filename)
        if d["atom_io"] is None:
            d["atom_io"] = {}
        return cls(**d)
