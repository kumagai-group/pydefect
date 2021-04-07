# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import argparse
import os
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.band_edge_states import BandEdgeStates
from pymatgen.core.structure import Structure
from vise.defaults import defaults
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.task import Task
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.util.file_transfer import FileLink


def is_file(filename):
    return Path(filename).is_file() and os.stat(filename).st_size != 0


def make_parchg_dir(band_edge_states: BandEdgeStates):

    if is_file("WAVECAR") is False:
        raise FileNotFoundError("WAVECAR does not exist or is empty.")

    Path("parchg").mkdir()
    os.chdir("parchg")
    user_incar_settings = {"LPARD": True,
                           "LSEPB": True,
                           "KPAR": 1,
                           "IBAND": band_edge_states.band_indices}

    options = CategorizedInputOptions(
        structure=Structure.from_file(Path.cwd() / ".." / defaults.contcar),
        task=Task.defect, xc=defaults.xc)

    vif = VaspInputFiles(options, user_incar_settings)
    vif.create_input_files(Path.cwd())
    FileLink(Path("../WAVECAR")).transfer(Path.cwd())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--band_edge_states",
                        default="band_edge_states.json", metavar="FILE",
                        type=loadfn)
    args = parser.parse_args()

    make_parchg_dir(args.band_edge_states)


if __name__ == "__main__":
    main()


