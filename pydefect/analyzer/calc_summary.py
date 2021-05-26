# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from dataclasses import dataclass
from typing import Optional, Dict

from monty.json import MSONable
from pydefect.analyzer.defect_structure_info import SymmRelation, DefectType
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class SingleCalcSummary(MSONable):
    charge: int
    atom_io: dict
    electronic_conv: Optional[bool] = None
    ionic_conv: Optional[bool] = None
    is_energy_strange: Optional[bool] = None
    same_config_from_init: Optional[bool] = None
    defect_type: Optional[str] = None
    symm_relation: Optional[str] = None
    donor_phs: Optional[bool] = None
    acceptor_phs: Optional[bool] = None
    unoccupied_deep_state: Optional[bool] = None
    occupied_deep_state: Optional[bool] = None
    same_structure: Optional[str] = None

    def same_atom_charge_io(self, other: "SingleCalcSummary"):
        return self.charge == other.charge and self.atom_io == other.atom_io

    @property
    def is_converged(self):
        return self.electronic_conv and self.ionic_conv

    @property
    def is_proper_result(self):
        return self.is_converged is True and self.is_energy_strange is False

    @property
    def is_unusual(self):
        return (self.symm_relation is SymmRelation.supergroup or
                self.same_config_from_init is False)

    @property
    def config_list(self):
        if self.is_proper_result:
            result = [".", ".", "."]
            if self.same_config_from_init in (False, None):
                result.extend([str(self.same_config_from_init),
                               str(self.defect_type),
                               str(self.symm_relation)])
            else:
                result.extend([".", ".", "."])
            return result
        elif self.is_converged:
            return [".", ".", self.is_energy_strange]
        else:
            return [self.electronic_conv, self.ionic_conv]


@dataclass
class CalcSummary(MSONable, ToJsonFileMixIn):
    single_summaries: Dict[str, SingleCalcSummary]

    def __str__(self):
        lines = [["name", "Ele. conv.", "Ionic conv.", "Is energy strange",
                  "Same config.", "Defect type", "Symm. Relation"]]
        for k, v in self.single_summaries.items():
            lines.append([k] + v.config_list)
        return tabulate(lines, stralign="center", tablefmt="pipe")
