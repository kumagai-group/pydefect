# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import fire
from monty.serialization import loadfn
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_structure_info import DefectStructureInfo, \
    fold_coords_in_structure
from pydefect.analyzer.vesta.create_vesta_file import VestaFile
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.supercell_info import SupercellInfo
from pymatgen.core import Structure, DummySpecies, Element


def make_defect_vesta_file(defect_structure: Structure,
                           def_str_info: DefectStructureInfo,
                           title: str = None,
                           cutoff: float = 5.0,
                           show_disp_min: float = 0.1,
                           disp_enhance_factor: float = 3.0,
                           initial_filename: str = "defect_initial.vesta",
                           final_filename: str = "defect.vesta"):
    species, initial_coords, final_coords, site_info = [], [], [], []
    vectors, vector_colors = {}, []
    i = 1
    for f_idx, (s, disp) in enumerate(zip(defect_structure, def_str_info.displacements)):
        if disp and disp.distance_from_defect < cutoff:
            species.append(disp.specie)
            initial_coords.append(disp.original_pos)
            final_coords.append(s.frac_coords)
            info = f"{round(disp.distance_from_defect, 1)}"

            if disp and disp.displace_distance > show_disp_min:
                vectors[i] = (v * disp_enhance_factor for v in disp.disp_vector)
                # use 200 as max instead of 255 since (255, 255, 255) is white.
                rgb = int(200 * disp.angle / 180)
                vector_colors.append((rgb, rgb, 100))
                info += f"_{round(disp.displace_distance, 2)}"

            site_info.append(info)
            i += 1

    initial_local = Structure(defect_structure.lattice, species, initial_coords)
    final_local = Structure(defect_structure.lattice, species, final_coords)
    for s, i in zip(final_local, site_info):
        s.properties["name"] = i
    for _, specie, fcoord in def_str_info.site_diff.inserted:
        final_local.append(Element(specie), coords=fcoord)
        final_local[-1].properties["name"] = f"{specie}_i"

    for (_, i_specie, i_fcoord), (_, r_specie, r_fcoord) in zip(def_str_info.site_diff.inserted_by_sub, def_str_info.site_diff.removed_by_sub):
        initial_local.append(Element(r_specie), coords=r_fcoord)
        final_local.append(Element(i_specie), coords=i_fcoord)
        final_local[-1].properties["name"] = f"{i_specie}_{r_specie}"

    for _, specie, fcoord in def_str_info.site_diff.removed:
        initial_local.append(Element(specie), coords=fcoord)
        final_local.append(DummySpecies(), coords=fcoord)
        final_local[-1].properties["name"] = f"Va_{specie}"

    if def_str_info.site_diff.is_complex_defect:
        final_local.append(DummySpecies(), coords=def_str_info.center)
        final_local[-1].properties["name"] = f"center"

    fold_coords_in_structure(initial_local, def_str_info.center)
    fold_coords_in_structure(final_local, def_str_info.center)

    boundary = sum([[x - 0.5, x + 0.5] for x in def_str_info.center], [])

    VestaFile(initial_local, title, boundary=boundary).write_file(initial_filename)
    VestaFile(final_local, title, vectors=vectors, vector_colors=vector_colors, boundary=boundary).write_file(final_filename)


def make_defect_vesta_file_via_fire(supercell_info_json,
                                    defect_entry_json,
                                    calc_results_json,
                                    defect_str_info_json,
                                    cutoff: float = 5.0,
                                    show_disp_min: float = 0.1,
                                    disp_enhance_factor: float = 3.0,
                                    initial_filename: str = "defect_initial.vesta",
                                    final_filename: str = "defect.vesta"):
    supercell_info: SupercellInfo = loadfn(supercell_info_json)
    defect_entry: DefectEntry = loadfn(defect_entry_json)
    calc_results: CalcResults = loadfn(calc_results_json)
    defect_str_info = loadfn(defect_str_info_json)
    title = f"{defect_entry.full_name} in {supercell_info.structure.formula}"
    make_defect_vesta_file(calc_results.structure,
                           defect_str_info,
                           title,
                           cutoff,
                           show_disp_min,
                           disp_enhance_factor,
                           initial_filename,
                           final_filename)


if __name__ == '__main__':
    fire.Fire(make_defect_vesta_file_via_fire)
