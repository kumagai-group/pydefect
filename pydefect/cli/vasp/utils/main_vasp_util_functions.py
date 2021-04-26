# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from pydefect.analyzer.grids import Grids
from pydefect.cli.vasp.make_defect_charge_info import make_spin_charges, \
    make_defect_charge_info
from pydefect.cli.vasp.make_light_chgcar import make_light_chgcar
from pymatgen.io.vasp import Chgcar
from vise.util.logger import get_logger

logger = get_logger(__name__)


def calc_defect_charge_info(args):
    band_idxs = []
    parchgs = []
    for parchg in args.parchgs:
        band_idx = int(parchg.split(".")[-2])
        logger.info(f"band index {band_idx} is being parsed.")
        # Use the band indices used in VASP, where it begins from 1.
        band_idxs.append(band_idx)
        p = Chgcar.from_file(parchg)
        parchgs.append(p)
        parent = Path(parchg).parent
        for spin, charge in zip(["up", "down"], make_spin_charges(p)):
            to_vesta_file = Path(f"defect_{band_idx}_{spin}.vesta")
            make_light_chgcar(charge,
                              parent / f"PARCHG_{band_idx}_{spin}",
                              vesta_file=args.vesta_file,
                              to_vesta_file=to_vesta_file)

    grids = None
    if args.grids_dirname:
        if (args.grids_dirname / "grids.npz").is_file():
            grids = Grids.from_file(args.grids_dirname / "grids.npz")
        else:
            grids = Grids.from_chgcar(parchgs[0])
            grids.dump(args.grids_dirname / "grids.npz")

    defect_charge_info = make_defect_charge_info(
        parchgs, band_idxs, args.bin_interval, grids)
    defect_charge_info.to_json_file()
    plt = defect_charge_info.show_dist()
    plt.savefig("dist.pdf")
