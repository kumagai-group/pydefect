# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, MutableMapping

from monty.json import MSONable
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pymatgen.symmetry.groups import SpaceGroup
from vise.util.mix_in import ToYamlFileMixIn
from vise.util.structure_symmetrizer import num_symmetry_operation


@dataclass
class Degeneracy(MSONable):
    site: int
    spin: int
    initial_site_sym: str = None
    final_site_sym: str = None

    @property
    def degeneracy(self) -> int:
        return self.site * self.spin


@dataclass
class Degeneracies(MutableMapping, ToYamlFileMixIn):
    d: Dict[str, Dict[int, Degeneracy]]

    def __iter__(self):
        return self.d.__iter__()

    def __len__(self) -> int:
        return len(self.d)

    def __getitem__(self, k):
        return self.d[k]

    def __delitem__(self, v) -> None:
        self.d.pop(v)

    def __setitem__(self, k, v) -> None:
        self.d[k] = v

    def as_dict(self):
        result = {}
        for k, v in self.d.items():
            result[k] = {}
            for kk, vv in v.items():
                vvv = vv.as_dict()
                vvv.pop("@class")
                vvv.pop("@module")
                vvv.pop("@version")
                result[k][kk] = vvv
        return result

    @classmethod
    def from_dict(cls, d):
        result = cls(d)
        for k, v in result.items():
            for kk, vv in v.items():
                result[k][kk] = Degeneracy.from_dict(vv)
        return result


class MakeDegeneracy:
    def __init__(self,
                 primitive_sg_symbol: str,
                 int_threshold: float = 0.1):
        self._int_threshold = int_threshold
        self._primitive_num_sym_opt = len(SpaceGroup(primitive_sg_symbol))
        self._deg_dict = defaultdict(dict)

    def add_degeneracy(self,
                       energy_info: DefectEnergyInfo,
                       calc_results: CalcResults,
                       structure_info: DefectStructureInfo):
        spin = self.mag_to_spin_degeneracy(calc_results.magnetization)
        site = (self._primitive_num_sym_opt
                / num_symmetry_operation(structure_info.final_site_sym))

        degeneracy = Degeneracy(int(site), spin,
                                structure_info.initial_site_sym,
                                structure_info.final_site_sym)
        self._deg_dict[energy_info.name][energy_info.charge] = degeneracy

    def mag_to_spin_degeneracy(self, mag: float) -> int:
        rounded_mag = round(mag)
        if abs(rounded_mag - mag) > self._int_threshold:
            raise ValueError
        return 2 * abs(rounded_mag) + 1

    @property
    def degeneracies(self):
        return Degeneracies(self._deg_dict)
