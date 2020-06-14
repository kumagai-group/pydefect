# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pymatgen import DummySpecie
from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
from pymatgen.io.vasp import Chgcar
from vise.util.structure_symmetrizer import StructureSymmetrizer


def interstitials_from_charge_density(
        chgcar: Chgcar,
        threshold_frac: float = None,
        threshold_abs: float = None,
        min_dist: float = 0.5,
        tol: float = 0.2,
        radius: float = 0.4) -> None:
    cda = ChargeDensityAnalyzer(chgcar=chgcar)
    cda.get_local_extrema(threshold_frac=threshold_frac,
                          threshold_abs=threshold_abs)
    cda.sort_sites_by_integrated_chg(r=radius)
    # Remove sites near host atoms.
    cda.remove_collisions(min_dist)
    # Cluster interstitials that are too close together using a tol.
    cda.cluster_nodes(tol=tol)

    structure = chgcar.structure.copy()
    print(cda.extrema_df)

    start_index = len(structure)
    end_index = len(structure) + len(cda.extrema_coords)
    interstitial_indices = [i for i in range(start_index, end_index)]
    coords = cda.extrema_coords
    for c in coords:
        structure.append(DummySpecie(), c)
    sga = StructureSymmetrizer(structure)
    equiv_atoms = sga.spglib_sym_data["equivalent_atoms"]

    print("")
    print("++ Inequivalent indices and site symmetries ++")
    orig_num_atoms = len(chgcar.structure)
    for i, ii in enumerate(interstitial_indices):
        if ii == equiv_atoms[ii]:
            idx = orig_num_atoms + i
            coords = sga.primitive[idx].frac_coords
            idx_coords = \
                f"{i:>3} {coords[0]:8.4f} {coords[1]:8.4f} {coords[2]:8.4f}"

            print(idx_coords, sga.spglib_sym_data["site_symmetry_symbols"][ii])
