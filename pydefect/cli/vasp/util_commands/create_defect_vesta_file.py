# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import fire
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.analyzer.vesta.create_vesta_file import VestaFile
from pymatgen import Structure, DummySpecies, loadfn


def make_defect_vesta_file(defect_structure: Structure,
                           def_str_info: DefectStructureInfo,
                           cutoff: float = 5.0,
                           show_disp_min: float = 0.1,
                           disp_enhance_factor: float = 3.0,
                           filename: str = "defect.vesta"):
    sites, vectors = [], {}
    i = 1
    for s, disp in zip(defect_structure, def_str_info.displacements):
        if disp and disp.distance_from_defect < cutoff:
            sites.append(s)
            if disp and disp.displace_distance > show_disp_min:
                vectors[i] = (v * disp_enhance_factor for v in disp.disp_vector)
            i += 1

    local = Structure.from_sites(sites)
    for _, _, fcoord in def_str_info.site_diff.removed:
        local.append(DummySpecies("X"), coords=fcoord)

    boundary = sum([[x - 0.5, x + 0.5] for x in def_str_info.center], [])
    vesta_file = VestaFile(local, vectors=vectors, boundary=boundary)
    vesta_file.write_file(filename)


def make_defect_vesta_file_via_fire(calc_results_json,
                                    defect_str_info_json,
                                    cutoff: float = 5.0,
                                    show_disp_min: float = 0.1,
                                    disp_enhance_factor: float = 3.0,
                                    filename: str = "defect.vesta"):
    calc_results = loadfn(calc_results_json)
    defect_str_info = loadfn(defect_str_info_json)
    make_defect_vesta_file(calc_results.structure, defect_str_info, cutoff,
                           show_disp_min, disp_enhance_factor, filename)


if __name__ == '__main__':
    fire.Fire(make_defect_vesta_file_via_fire)
