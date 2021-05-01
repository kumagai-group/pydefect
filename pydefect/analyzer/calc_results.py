# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from dataclasses import dataclass
from typing import List, Optional

from monty.json import MSONable
from pymatgen.core import IStructure
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn


logger = get_logger(__name__)


@dataclass
class CalcResults(MSONable, ToJsonFileMixIn):
    # keep structure and site_symmetry in CalcResults
    structure: IStructure
    energy: float
    magnetization: float
    # potential acting on the positive unit charge, whose sign is reserved from
    # vasp convention of atomic site potential.
    potentials: List[float]
    electronic_conv: Optional[bool] = None
    ionic_conv: Optional[bool] = None

    def show_convergence_warning(self):
        if self.electronic_conv is False:
            logger.warning("SCF is not converged.")
        if self.ionic_conv is False:
            logger.warning("Ionic relaxation is not converged.")

    def __str__(self):
        return f""" -- calc results info
energy: {self.energy:10.3f}
magnetization: {self.magnetization:6.2f}
electronic convergence: {self.electronic_conv}
ionic convergence: {self.ionic_conv}"""


class NoElectronicConvError(AssertionError):
    pass


class NoIonicConvError(AssertionError):
    pass
