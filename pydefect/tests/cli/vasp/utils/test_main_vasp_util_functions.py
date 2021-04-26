# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.calc_results import CalcResults
from pydefect.chem_pot_diag.chem_pot_diag import CompositionEnergies, \
    CompositionEnergy
from pydefect.cli.main_functions import calc_defect_structure_info
from pydefect.corrections.efnv_correction import ExtendedFnvCorrection
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.supercell_info import SupercellInfo
from pymatgen.core import Composition


