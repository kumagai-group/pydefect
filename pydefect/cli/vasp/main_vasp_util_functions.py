# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import os
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.grids import Grids
from pydefect.analyzer.refine_defect_structure import refine_defect_structure
from pydefect.cli.vasp.make_defect_charge_info import make_defect_charge_info
from pydefect.cli.vasp.get_defect_charge_state import get_defect_charge_state
from pydefect.input_maker.defect_entry import make_defect_entry
from pymatgen.core import Structure
from pymatgen.io.vasp import Chgcar
from vise.input_set.incar import ViseIncar
from vise.util.file_transfer import FileLink
from vise.util.logger import get_logger
from pymatgen.io.vasp.inputs import Poscar, Incar, Potcar

logger = get_logger(__name__)


def is_file(filename):
    return Path(filename).is_file() and os.stat(filename).st_size != 0


def calc_charge_state(args):
    poscar = Poscar.from_file(args.dir / "POSCAR")
    potcar = Potcar.from_file(args.dir / "POTCAR")
    incar = Incar.from_file(args.dir / "INCAR")
    charge_state = get_defect_charge_state(poscar, potcar, incar)
    logger.info(f"Charge state in {args.dir} is {charge_state}.")
    return charge_state


def make_defect_entry_main(args):
    charge_state = calc_charge_state(args)
    structure = Structure.from_file(args.dir / "POSCAR")
    defect_entry = make_defect_entry(name=args.name,
                                     charge=charge_state,
                                     perfect_structure=args.perfect,
                                     defect_structure=structure)
    defect_entry.to_json_file()


def make_parchg_dir(args):
    os.chdir(args.dir)
    if is_file("WAVECAR") is False:
        raise FileNotFoundError("WAVECAR does not exist or is empty.")

    try:
        calc_results: CalcResults = loadfn("calc_results.json")
    except FileNotFoundError:
        logger.info("Need to create calc_results.json beforehand.")
        raise
    calc_results.show_convergence_warning()

    # Increment index by 1 as VASP band index begins from 1.
    incar = ViseIncar.from_file("INCAR")
    band_edge_states = loadfn("band_edge_states.json")
    iband = [i + 1 for i in band_edge_states.band_indices_from_vbm_to_cbm]
    incar.update({"LPARD": True, "LSEPB": True, "KPAR": 1, "IBAND": iband})

    parchg = Path("parchg")
    parchg.mkdir()
    os.chdir("parchg")
    incar.write_file("INCAR")
    FileLink(Path("../WAVECAR")).transfer(Path.cwd())
    FileLink(Path("../POSCAR")).transfer(Path.cwd())
    FileLink(Path("../POTCAR")).transfer(Path.cwd())
    FileLink(Path("../KPOINTS")).transfer(Path.cwd())
    os.chdir("..")


def make_refine_defect_poscar(args):
    structure = refine_defect_structure(args.structure,
                                        args.defect_entry.anchor_atom_index,
                                        args.defect_entry.anchor_atom_coords)
    if structure:
        print(structure.to(fmt="poscar", filename=args.poscar_name))


def calc_grids(args):
    grids = Grids.from_chgcar(args.chgcar)
    grids.dump()


def make_defect_charge_info_main(args):
    band_idxs = [int(parchg.split(".")[-2]) - 1 for parchg in args.parchgs]
    parchgs = [Chgcar.from_file(parchg) for parchg in args.parchgs]
    defect_charge_info = make_defect_charge_info(
        parchgs, band_idxs, args.bin_interval, args.grids)
    defect_charge_info.to_json_file()
    plt = defect_charge_info.show_dist()
    plt.savefig("dist.pdf")


