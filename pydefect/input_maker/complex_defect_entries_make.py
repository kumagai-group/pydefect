# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
from pydefect.input_maker.complex_defect_set import ComplexDefectSet, \
    ComplexDefect
from pydefect.input_maker.defect_entries_maker import copy_to_structure, \
    to_istructure
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.supercell_info import SupercellInfo
from pymatgen.core import Element


class ComplexDefectEntriesMaker:
    def __init__(self,
                 supercell_info: SupercellInfo,
                 complex_defect_set: ComplexDefectSet):
        self.supercell_info = supercell_info
        self.defect_entries = list()

        for defect in complex_defect_set:
            (structure, coords, site_sym,
             perturbed_structure, perturbed_sites, perturbed_site_sym) = \
                self._create_defect_structures(defect)

            for charge in defect.charges:
                self.defect_entries.append(
                    DefectEntry(name=defect.name,
                                charge=charge,
                                structure=structure,
                                site_symmetry=site_sym,
                                defect_center=coords,
                                perturbed_structure=perturbed_structure,
                                perturbed_sites=perturbed_sites,
                                perturbed_site_symmetry=perturbed_site_sym))

    def _create_defect_structures(self, defect: ComplexDefect):
        structure = copy_to_structure(self.supercell_info.structure)
        for out_idx, in_elem in sorted(defect.replaced_sites.items(),
                                       reverse=True):
            if in_elem is None:
                structure.pop(out_idx)
            else:
                structure.sites[out_idx].species = Element(in_elem)

        coords = [0.125, 0.125, 0.125]
        site_symmetry = "m-3m"
        return (to_istructure(structure), tuple(coords), site_symmetry,
                None, None, None)
