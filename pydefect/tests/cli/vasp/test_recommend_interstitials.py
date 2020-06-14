# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from pydefect.cli.vasp.recommend_interstitials import \
    interstitials_from_charge_density


def test_interstitials_from_charge_density():
    interstitials_from_charge_density(str(Path(__file__).parent / "vasp_files" / "H2_CHGCAR"))


"""
TODO
-

DONE
"""