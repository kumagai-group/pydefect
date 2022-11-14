# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import itertools
from itertools import groupby
from typing import List

import numpy as np
import pandas as pd
from pandas import DataFrame
from pydefect.analyzer.defect_structure_info import remove_dot
from pydefect.input_maker.local_extrema import VolumetricDataLocalExtrema, \
    CoordInfo, VolumetricDataAnalyzeParams
from pydefect.util.structure_tools import Distances
from pymatgen.core import Element, Structure
from pymatgen.io.vasp import VolumetricData, Chgcar
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
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


class ChargeDensityAnalyzer:
    """
    Analyzer to find potential interstitial sites based on charge density. The
    `total` charge density is used.
    """

    def __init__(self, chgcar):
        """
        Initialization.

        Args:
            chgcar (pmg.Chgcar): input Chgcar object.
        """
        self.chgcar = chgcar
        self.structure = chgcar.structure
        self.extrema_coords = []  # list of frac_coords of local extrema
        self.extrema_type = None  # "local maxima" or "local minima"
        self._extrema_df = None  # extrema frac_coords - chg density table
        self._charge_distribution_df = None  # frac_coords - chg density table

    @classmethod
    def from_file(cls, chgcar_filename):
        """
        Init from a CHGCAR.

        :param chgcar_filename:
        :return:
        """
        chgcar = Chgcar.from_file(chgcar_filename)
        return cls(chgcar=chgcar)

    @property
    def charge_distribution_df(self):
        """
        :return: Charge distribution.
        """
        if self._charge_distribution_df is None:
            return self._get_charge_distribution_df()
        return self._charge_distribution_df

    @property
    def extrema_df(self):
        """
        :return: The extrema in charge density.
        """
        if self.extrema_type is None:
            logger.warning("Please run ChargeDensityAnalyzer.get_local_extrema first!")
        return self._extrema_df

    def _get_charge_distribution_df(self):
        """
        Return a complete table of fractional coordinates - charge density.
        """
        # Fraction coordinates and corresponding indices
        axis_grid = np.array([np.array(self.chgcar.get_axis_grid(i)) / self.structure.lattice.abc[i] for i in range(3)])
        axis_index = np.array([range(len(axis_grid[i])) for i in range(3)])

        data = {}

        for index in itertools.product(*axis_index):
            a, b, c = index
            f_coords = (axis_grid[0][a], axis_grid[1][b], axis_grid[2][c])
            data[f_coords] = self.chgcar.data["total"][a][b][c]

        # Fraction coordinates - charge density table
        df = pd.Series(data).reset_index()
        df.columns = ["a", "b", "c", "Charge Density"]
        self._charge_distribution_df = df

        return df

    def _update_extrema(self, f_coords, extrema_type, threshold_frac=None, threshold_abs=None):
        """Update _extrema_df, extrema_type and extrema_coords"""

        if threshold_frac is not None:
            if threshold_abs is not None:
                logger.warning("Filter can be either threshold_frac or threshold_abs!")  # Exit if both filter are set
                return
            if threshold_frac > 1 or threshold_frac < 0:
                raise Exception("threshold_frac range is [0, 1]!")

        # Return empty result if coords list is empty
        if len(f_coords) == 0:
            df = pd.DataFrame({}, columns=["A", "B", "C", "Chgcar"])
            self._extrema_df = df
            self.extrema_coords = []
            logger.info(f"Find {len(df)} {extrema_type}.")
            return

        data = {}
        unit = 1 / np.array(self.chgcar.dim)  # pixel along a, b, c

        for fc in f_coords:
            a, b, c = tuple(map(int, fc / unit))
            data[tuple(fc)] = self.chgcar.data["total"][a][b][c]

        df = pd.Series(data).reset_index()
        df.columns = ["a", "b", "c", "Charge Density"]
        ascending = extrema_type == "local minima"

        if threshold_abs is None:
            threshold_frac = threshold_frac if threshold_frac is not None else 1.0
            num_extrema = int(threshold_frac * len(f_coords))
            df = df.sort_values(by="Charge Density", ascending=ascending)[0:num_extrema]
            df.reset_index(drop=True, inplace=True)  # reset major index
        else:  # threshold_abs is set
            df = df.sort_values(by="Charge Density", ascending=ascending)
            df = df[df["Charge Density"] <= threshold_abs] if ascending else df[df["Charge Density"] >= threshold_abs]

        extrema_coords = []
        for row in df.iterrows():
            fc = np.array(row[1]["a":"c"])
            extrema_coords.append(fc)

        self._extrema_df = df
        self.extrema_type = extrema_type
        self.extrema_coords = extrema_coords
        logger.info(f"Find {len(df)} {extrema_type}.")

    def get_local_extrema(self, find_min=True, threshold_frac=None, threshold_abs=None):
        """
        Get all local extrema fractional coordinates in charge density,
        searching for local minimum by default. Note that sites are NOT grouped
        symmetrically.

        Args:
            find_min (bool): True to find local minimum else maximum, otherwise
                find local maximum.

            threshold_frac (float): optional fraction of extrema shown, which
                returns `threshold_frac * tot_num_extrema` extrema fractional
                coordinates based on highest/lowest intensity.

                E.g. set 0.2 to show the extrema with 20% highest or lowest
                intensity. Value range: 0 <= threshold_frac <= 1

                Note that threshold_abs and threshold_frac should not set in the
                same time.

            threshold_abs (float): optional filter. When searching for local
                minima, intensity <= threshold_abs returns; when searching for
                local maxima, intensity >= threshold_abs returns.

                Note that threshold_abs and threshold_frac should not set in the
                same time.

        Returns:
            extrema_coords (list): list of fractional coordinates corresponding
                to local extrema.
        """
        try:
            from skimage.feature import peak_local_max
        except ImportError:
            logger.warning("To find the extrema of the coordinates, "
                           "install skimage.")
            raise

        sign, extrema_type = 1, "local maxima"

        if find_min:
            sign, extrema_type = -1, "local minima"

        # Make 3x3x3 supercell
        # This is a trick to resolve the periodical boundary issue.
        total_chg = sign * self.chgcar.data["total"]
        total_chg = np.tile(total_chg, reps=(3, 3, 3))
        coordinates = peak_local_max(total_chg, min_distance=1)

        # Remove duplicated sites introduced by supercell.
        f_coords = [coord / total_chg.shape * 3 for coord in coordinates]
        f_coords = [f - 1 for f in f_coords if all(np.array(f) < 2) and all(np.array(f) >= 1)]

        # Update information
        self._update_extrema(
            f_coords,
            extrema_type,
            threshold_frac=threshold_frac,
            threshold_abs=threshold_abs,
        )

        return self.extrema_coords

    def cluster_nodes(self, tol=0.2):
        """
        Cluster nodes that are too close together using a tol.

        Args:
            tol (float): A distance tolerance. PBC is taken into account.
        """
        lattice = self.structure.lattice
        vf_coords = self.extrema_coords

        if len(vf_coords) == 0:
            if self.extrema_type is None:
                logger.warning("Please run ChargeDensityAnalyzer.get_local_extrema first!")
                return None
            new_f_coords = []
            self._update_extrema(new_f_coords, self.extrema_type)
            return new_f_coords

        # Manually generate the distance matrix (which needs to take into
        # account PBC.
        dist_matrix = np.array(lattice.get_all_distances(vf_coords, vf_coords))
        dist_matrix = (dist_matrix + dist_matrix.T) / 2

        for i in range(len(dist_matrix)):
            dist_matrix[i, i] = 0
        condensed_m = squareform(dist_matrix)
        z = linkage(condensed_m)
        cn = fcluster(z, tol, criterion="distance")
        merged_fcoords = []

        for n in set(cn):
            frac_coords = []
            for i, j in enumerate(np.where(cn == n)[0]):
                if i == 0:
                    frac_coords.append(self.extrema_coords[j])
                else:
                    f_coords = self.extrema_coords[j]
                    # We need the image to combine the frac_coords properly.
                    d, image = lattice.get_distance_and_image(frac_coords[0], f_coords)
                    frac_coords.append(f_coords + image)
            merged_fcoords.append(np.average(frac_coords, axis=0))

        merged_fcoords = [f - np.floor(f) for f in merged_fcoords]
        merged_fcoords = [f * (np.abs(f - 1) > 1e-15) for f in merged_fcoords]
        # the second line for fringe cases like
        # np.array([ 5.0000000e-01 -4.4408921e-17  5.0000000e-01])
        # where the shift to [0,1) does not work due to float precision
        self._update_extrema(merged_fcoords, extrema_type=self.extrema_type)
        logger.debug(f"{len(self.extrema_coords)} vertices after combination.")
        return None

    def remove_collisions(self, min_dist=0.5):
        """
        Remove predicted sites that are too close to existing atoms in the
        structure.

        Args:
            min_dist (float): The minimum distance (in Angstrom) that
                a predicted site needs to be from existing atoms. A min_dist
                with value <= 0 returns all sites without distance checking.
        """
        s_f_coords = self.structure.frac_coords
        f_coords = self.extrema_coords
        if len(f_coords) == 0:
            if self.extrema_type is None:
                logger.warning("Please run ChargeDensityAnalyzer.get_local_extrema first!")
                return None
            new_f_coords = []
            self._update_extrema(new_f_coords, self.extrema_type)
            return new_f_coords

        dist_matrix = self.structure.lattice.get_all_distances(f_coords, s_f_coords)
        all_dist = np.min(dist_matrix, axis=1)
        new_f_coords = []

        for i, f in enumerate(f_coords):
            if all_dist[i] > min_dist:
                new_f_coords.append(f)
        self._update_extrema(new_f_coords, self.extrema_type)

        return new_f_coords

    def get_structure_with_nodes(
            self,
            find_min=True,
            min_dist=0.5,
            tol=0.2,
            threshold_frac=None,
            threshold_abs=None,
    ):
        """
        Get the modified structure with the possible interstitial sites added.
        The species is set as a DummySpecies X.

        Args:
            find_min (bool): True to find local minimum else maximum, otherwise
                find local maximum.

            min_dist (float): The minimum distance (in Angstrom) that
                a predicted site needs to be from existing atoms. A min_dist
                with value <= 0 returns all sites without distance checking.

            tol (float): A distance tolerance of nodes clustering that sites too
                closed to other predicted sites will be merged. PBC is taken
                into account.

            threshold_frac (float): optional fraction of extrema, which returns
                `threshold_frac * tot_num_extrema` extrema fractional
                coordinates based on highest/lowest intensity.

                E.g. set 0.2 to insert DummySpecies atom at the extrema with 20%
                highest or lowest intensity.
                Value range: 0 <= threshold_frac <= 1

                Note that threshold_abs and threshold_frac should not set in the
                same time.

            threshold_abs (float): optional filter. When searching for local
                minima, intensity <= threshold_abs returns; when searching for
                local maxima, intensity >= threshold_abs returns.

                Note that threshold_abs and threshold_frac should not set in the
                same time.

        Returns:
            structure (Structure)
        """

        structure = self.structure.copy()
        self.get_local_extrema(
            find_min=find_min,
            threshold_frac=threshold_frac,
            threshold_abs=threshold_abs,
        )

        self.remove_collisions(min_dist)
        self.cluster_nodes(tol=tol)
        for fc in self.extrema_coords:
            structure.append("X", fc)

        return structure

    def sort_sites_by_integrated_chg(self, r=0.4):
        """
        Get the average charge density around each local minima in the charge density
        and store the result in _extrema_df
        Args:
            r (float): radius of sphere around each site to evaluate the average
        """

        if self.extrema_type is None:
            self.get_local_extrema()
        int_den = []
        for isite in self.extrema_coords:
            mask = self._dist_mat(isite) < r
            vol_sphere = self.chgcar.structure.volume * (mask.sum() / self.chgcar.ngridpts)
            chg_in_sphere = np.sum(self.chgcar.data["total"] * mask) / mask.size / vol_sphere
            int_den.append(chg_in_sphere)
        self._extrema_df["avg_charge_den"] = int_den
        self._extrema_df.sort_values(by=["avg_charge_den"], inplace=True)
        self._extrema_df.reset_index(drop=True, inplace=True)

    def _dist_mat(self, pos_frac):
        # return a matrix that contains the distances
        aa = np.linspace(0, 1, len(self.chgcar.get_axis_grid(0)), endpoint=False)
        bb = np.linspace(0, 1, len(self.chgcar.get_axis_grid(1)), endpoint=False)
        cc = np.linspace(0, 1, len(self.chgcar.get_axis_grid(2)), endpoint=False)
        AA, BB, CC = np.meshgrid(aa, bb, cc, indexing="ij")
        dist_from_pos = self.chgcar.structure.lattice.get_all_distances(
            fcoords1=np.vstack([AA.flatten(), BB.flatten(), CC.flatten()]).T,
            fcoords2=pos_frac,
        )
        return dist_from_pos.reshape(AA.shape)


def extrema_coords(volumetric_data: VolumetricData,
                   find_min: bool,
                   params: VolumetricDataAnalyzeParams) -> DataFrame:
    # try:
    #     from pymatgen.analysis.defects.utils import ChargeDensityAnalyzer
    # except ImportError:
    #     logger.warning("To find the extrema of the coordinates, "
    #                    "install pymatgen-analysis-defects code.")
    #     raise

    logger.info("When using this code, cite the paper related to skimage. "
                "see https://peerj.com/articles/453")

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

