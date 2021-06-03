# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from itertools import groupby
from typing import List

from pandas import DataFrame
from pydefect.analyzer.defect_structure_info import remove_dot
from pydefect.input_maker.local_extrema import VolumetricDataLocalExtrema, \
    CoordInfo, VolumetricDataAnalyzeParams
from pydefect.util.structure_tools import Distances
from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
from pymatgen.core import Element, Structure
from pymatgen.io.vasp import VolumetricData
from vise.util.logger import get_logger
from vise.util.structure_symmetrizer import StructureSymmetrizer

logger = get_logger(__name__)


def make_local_extrema_from_volumetric_data(
        volumetric_data: VolumetricData,
        params: VolumetricDataAnalyzeParams,
        info: str = None,
        find_min: bool = True) -> VolumetricDataLocalExtrema:
    extrema = extrema_coords(volumetric_data, find_min, params)
    coord_infos = find_inequivalent_coords(volumetric_data.structure,
                                           extrema)
    return VolumetricDataLocalExtrema(volumetric_data.structure,
                                      is_min=find_min,
                                      extrema_points=coord_infos,
                                      info=info,
                                      params=params)


def find_inequivalent_coords(structure: Structure,
                             df: DataFrame) -> List[CoordInfo]:
    result = []
    initial_sg = StructureSymmetrizer(structure).sg_number
    added_structure = Structure.from_dict(structure.as_dict())
    start_index = len(structure)
    for _, column in df.iterrows():
        coords = [column.a, column.b, column.c]
        assert Element.Og not in structure.composition.elements
        added_structure.append(Element.Og, coords)
    end_index = len(added_structure)

    symmetrizer = StructureSymmetrizer(added_structure)

    if initial_sg != symmetrizer.sg_number:
        logger.warning("The symmetry has changed, meaning all the symmetry "
                       "equivalent sites do not exist.")

    sym_data = symmetrizer.spglib_sym_data
    _indices = [i for i in range(start_index, end_index)]
    repr_atom_pairs = zip(sym_data["equivalent_atoms"][_indices], _indices)

    key = lambda x: x[0]
    for _, equiv_sites in groupby(sorted(repr_atom_pairs, key=key), key=key):
        coords, quantity = [], []
        for repr_idx, atom_idx in equiv_sites:
            fcoord = added_structure[atom_idx].frac_coords
            coords.append(tuple(fcoord))
            key = "ave_value" if "ave_value" in df else "value"
            quantity.append(df[key][atom_idx - start_index])
            if repr_idx == atom_idx:
                site_sym = sym_data["site_symmetry_symbols"][atom_idx]
                distances = Distances(structure, fcoord)
                coordination = distances.coordination()

        coord_info = CoordInfo(site_symmetry=remove_dot(site_sym),
                               coordination=coordination,
                               frac_coords=coords,
                               quantities=quantity)
        result.append(coord_info)
    return result


def extrema_coords(volumetric_data: VolumetricData,
                   find_min: bool,
                   params: VolumetricDataAnalyzeParams) -> DataFrame:
    result = ChargeDensityAnalyzer(chgcar=volumetric_data)
    result.get_local_extrema(threshold_frac=params.threshold_frac,
                             threshold_abs=params.threshold_abs,
                             find_min=find_min)
    if params.min_dist:
        # Remove sites near host atoms.
        result.remove_collisions(params.min_dist)
    if params.tol:
        try:
            # Cluster interstitials that are too close together using a tol.
            result.cluster_nodes(tol=params.tol)
        except ValueError:
            # error raised when only a single local extreme point.
            pass
    if params.radius:
        result.sort_sites_by_integrated_chg(r=params.radius)
    result = result.extrema_df
    if "avg_charge_den" in result.columns:
        result.rename(columns={'avg_charge_den': 'ave_value'}, inplace=True)
    result.rename(columns={'Charge Density': 'value'}, inplace=True)
    logger.info("\n" + result.__str__())
    return result

