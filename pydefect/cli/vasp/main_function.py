# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pydefect.input_maker.supercell_maker import SupercellMaker


def make_supercell(args):
    maker = SupercellMaker(args.unitcell, args.matrix)
    maker.supercell.structure.to(filename="SPOSCAR")
