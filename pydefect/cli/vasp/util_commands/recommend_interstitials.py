# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.

import fire
from pymatgen import Element, Structure
from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
from pymatgen.io.vasp import Chgcar, VolumetricData
from vise.util.structure_symmetrizer import StructureSymmetrizer


def interstitials_from_volumetric_data(
        volumetric_data: VolumetricData,
        find_min: bool = True,
        threshold_frac: float = None,
        threshold_abs: float = None,
        min_dist: float = 0.5,
        tol: float = 0.5,
        radius: float = 0.4) -> None:

    cda = ChargeDensityAnalyzer(chgcar=volumetric_data)
    cda.get_local_extrema(threshold_frac=threshold_frac,
                          threshold_abs=threshold_abs,
                          find_min=find_min)
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

    structure: Structure = volumetric_data.structure.copy()
    start_index = len(structure)
    end_index = len(structure) + len(cda.extrema_coords)
    interstitial_indices = [i for i in range(start_index, end_index)]
    for coords in cda.extrema_coords:
        structure.append(Element.H, coords)
    symmetrizer = StructureSymmetrizer(structure)
    equiv_atoms = symmetrizer.spglib_sym_data["equivalent_atoms"]

    print(f"Host symmetry {symmetrizer.spglib_sym_data['international']}")
    print("++ Inequivalent indices and site symmetries ++")
    orig_num_atoms = len(volumetric_data.structure)

    structure_for_visualize: Structure = volumetric_data.structure.copy()

    for i, ii in enumerate(interstitial_indices):
        if ii == equiv_atoms[ii]:
            idx = orig_num_atoms + i
            coords = structure[idx].frac_coords

            def element():
                for z in range(1, 100):
                    e = Element.from_Z(z)
                    if e not in structure_for_visualize.composition:
                        yield e

            structure_for_visualize.append(next(element()), coords)

            idx_coords = \
                f"{i:>3} {coords[0]:8.4f} {coords[1]:8.4f} {coords[2]:8.4f}"
            print(idx_coords,
                  symmetrizer.spglib_sym_data["site_symmetry_symbols"][ii])

        structure_for_visualize.to(filename="POSCAR-w-interstitials")


def interstitials_from_charge_density(aeccar0, aeccar2, **kwargs):
    aeccar = Chgcar.from_file(aeccar0) + Chgcar.from_file(aeccar2)
    interstitials_from_volumetric_data(aeccar, **kwargs)


if __name__ == '__main__':
    fire.Fire(interstitials_from_charge_density)
