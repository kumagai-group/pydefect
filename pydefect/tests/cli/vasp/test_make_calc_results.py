# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
from pymatgen import Structure
from pymatgen.io.vasp import Vasprun, Outcar

from pydefect.cli.vasp.make_calc_results import MakeCalcResultsFromVaspResults


def test_make_calc_results_from_vasp_results(vasp_files):
    vasprun = Vasprun(vasp_files / "MgO_conv_Va_O_0" / "vasprun.xml")
    outcar = Outcar(vasp_files / "MgO_conv_Va_O_0" / "OUTCAR")
    results = MakeCalcResultsFromVaspResults(vasprun, outcar)

    expected_structure = \
        Structure.from_file(vasp_files / "MgO_conv_Va_O_0" / "CONTCAR")
    assert results.calc_results.structure == expected_structure
    assert results.calc_results.site_symmetry == "m-3m"
    assert results.calc_results.energy == -34.91084360
    assert results.calc_results.magnetization == 1.03e-05

"""
TODO
-
        # vasprun related
        vasprun = parse_file(Vasprun, Path(directory_path) / vasprun)
        if vasprun.converged_electronic is False:
            raise NoConvergenceError("Electronic step is not converged.")
        if vasprun.converged_ionic is False:
            logger.warning("Ionic step is not converged.")

        # outcar related
        outcar = parse_file(Outcar, Path(directory_path) / outcar)
        magnetization = outcar.total_mag or 0.0

        try:
            _, vbm, cbm = \
                (i["energy"] for i in band_gap_properties(vasprun, outcar))
        except TypeError:
            vbm = cbm = vasprun.efermi

        # contcar related
        contcar = parse_file(Poscar.from_file, Path(directory_path) / contcar)
        final_structure = contcar.structure
        volume = contcar.structure.volume
        sga = SpacegroupAnalyzer(
            final_structure, defect_symprec, angle_tolerance)
        site_symmetry = sga.get_point_group_symbol()

        # If defect_entry is None, system is regarded as perfect supercell.
        if not defect_entry:
            center = neighboring_sites = defect_coords = displacements = None
        else:
            if defect_entry.defect_type.is_defect_center_atom:
                center = defect_entry.inserted_atoms[0]["index"]
                defect_coords = list(final_structure[center].frac_coords)
            else:
                center = defect_coords = defect_entry.defect_center_coords

            if cutoff is None:
                min_d = min_distance_from_coords(final_structure, defect_coords)
                cutoff = round(min_d * CUTOFF_FACTOR, 2)

            neighboring_sites = []
            for i, site in enumerate(final_structure):
                # Calculate the distance between defect and site
                d, _ = site.distance_and_image_from_frac_coords(defect_coords)
                # Defect itself is excluded from neighboring sites in Pydefect.
                if 1e-5 < d < cutoff:
                    neighboring_sites.append(i)
            if not neighboring_sites:
                raise StructureError(f"No neighboring site detected. Increase "
                                     f"cutoff radius from {cutoff}.")

            displacements = get_displacements(final_structure,
                                              defect_entry.initial_structure,
                                              center,
                                              defect_entry.anchor_atom_index)
        # procar related
        kwargs = {}
        if procar:
            procar = parse_file(Procar, Path(directory_path) / procar)
            # The k-point indices at the band edges in defect calculations.
            # hob (lub) = highest (lowest) (un)occupied state
            # The small number (0.1) must be added to avoid magnetization = 0.0
            up_index = round((outcar.nelect + (magnetization + 0.1)) / 2) - 1
            down_index = round((outcar.nelect - (magnetization + 0.1)) / 2) - 1
            hob_index = {Spin.up: up_index, Spin.down: down_index}

            prop = ProcarDefectProperty. \
                analyze_procar(hob_index, procar, vasprun.eigenvalues,
                               final_structure, neighboring_sites)
            for k, v in prop.as_dict().items():
                if k[0] != "@":
                    kwargs[k] = v
DONE
"""