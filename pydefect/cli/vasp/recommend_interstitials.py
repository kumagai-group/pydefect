# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import fire
from pymatgen import Element
from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
from pymatgen.io.vasp import Chgcar
from vise.util.structure_symmetrizer import StructureSymmetrizer


def interstitials_from_charge_density(
        chgcar_file: str,
        threshold_frac: float = None,
        threshold_abs: float = None,
        min_dist: float = 0.5,
        tol: float = 0.2,
        radius: float = 0.4) -> None:
    chgcar = Chgcar.from_file(chgcar_file)
    cda = ChargeDensityAnalyzer(chgcar=chgcar)
    cda.get_local_extrema(threshold_frac=threshold_frac,
                          threshold_abs=threshold_abs)
    cda.sort_sites_by_integrated_chg(r=radius)
    # Remove sites near host atoms.
    cda.remove_collisions(min_dist)
    # Cluster interstitials that are too close together using a tol.
    try:
        cda.cluster_nodes(tol=tol)
    except ValueError:
        # error raised when only a single local extreme point.
        pass
    print(cda.extrema_df)

    structure = chgcar.structure.copy()
    start_index = len(structure)
    end_index = len(structure) + len(cda.extrema_coords)
    interstitial_indices = [i for i in range(start_index, end_index)]
    for coords in cda.extrema_coords:
        structure.append(Element.H, coords)
    symmetrizer = StructureSymmetrizer(structure)
    equiv_atoms = symmetrizer.spglib_sym_data["equivalent_atoms"]

    print("")
    print("++ Inequivalent indices and site symmetries ++")
    orig_num_atoms = len(chgcar.structure)
    for i, ii in enumerate(interstitial_indices):
        if ii == equiv_atoms[ii]:
            idx = orig_num_atoms + i
            coords = symmetrizer.primitive[idx].frac_coords
            idx_coords = \
                f"{i:>3} {coords[0]:8.4f} {coords[1]:8.4f} {coords[2]:8.4f}"

            print(idx_coords, symmetrizer.spglib_sym_data["site_symmetry_symbols"][ii])


if __name__ == '__main__':
    fire.Fire(interstitials_from_charge_density)
