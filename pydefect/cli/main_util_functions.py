# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.cli.main_functions import get_calc_results, parse_dirs
from pydefect.cli.make_defect_vesta_file import MakeDefectVestaFile
from pydefect.cli.vasp.make_composition_energies_from_mp import \
    make_composition_energies_from_mp
from pydefect.cli.vasp.make_gkfo_correction import make_gkfo_correction
from pydefect.corrections.site_potential_plotter import SitePotentialMplPlotter
from vise.util.logger import get_logger


logger = get_logger(__name__)


def composition_energies_from_mp(args) -> None:
    composition_energies = make_composition_energies_from_mp(
        elements=args.elements, atom_energy_yaml=args.atom_energy_yaml)
    composition_energies.to_yaml_file()


def add_interstitials_from_local_extrema(args) -> None:
    supercell_info = args.local_extrema.append_sites_to_supercell_info(
        args.supercell_info, args.indices)
    supercell_info.to_json_file()


def make_defect_vesta_file(args) -> None:
    def _inner(_dir: Path):
        defect_str_info = loadfn(_dir / "defect_structure_info.json")
        make_vesta_file = MakeDefectVestaFile(defect_str_info,
                                              args.cutoff,
                                              args.min_displace_w_arrows,
                                              args.arrow_factor,
                                              args.title)
        make_vesta_file.initial_vesta.write_file(_dir / "defect_initial.vesta")
        make_vesta_file.final_vesta.write_file(_dir / "defect.vesta")

    parse_dirs(args.dirs, _inner)


def make_gkfo_correction_from_vasp(args):
    gkfo = make_gkfo_correction(
        efnv_correction=args.initial_efnv_correction,
        additional_charge=args.charge_diff,
        final_calc_results=args.final_calc_results,
        initial_calc_results=args.initial_calc_results,
        diele_tensor=args.unitcell.dielectric_constant,
        ion_clamped_diele_tensor=args.unitcell.ele_dielectric_const)
    print(gkfo)
    gkfo.to_json_file("gkfo_correction.json")
    plotter = SitePotentialMplPlotter.from_gkfo_corr(
        title="GKFO correction", gkfo_correction=gkfo)
    plotter.construct_plot()
    plotter.plt.savefig(fname="gkfo_correction.pdf")
    plotter.plt.clf()


