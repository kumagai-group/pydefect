# -*- coding: utf-8 -*-
#  Copyright (c) 2023 Kumagai group.
import pytest
from pydefect.input_maker.complex_defect_entries_make import \
    ComplexDefectEntriesMaker
from pydefect.input_maker.complex_defect_set import ComplexDefectSet, \
    ComplexDefect
from pydefect.input_maker.defect_entry import DefectEntry
from pydefect.input_maker.supercell_info import Interstitial, SupercellInfo
from pymatgen.core import Lattice, IStructure
from vise.util.structure_symmetrizer import Site


@pytest.fixture
def supercell_info():
    sites = {"H1": Site(element="H", wyckoff_letter="a", site_symmetry="m-3m",
                        equivalent_atoms=[0]),
             "He1": Site(element="He", wyckoff_letter="b", site_symmetry="m-3m",
                         equivalent_atoms=[1])}
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]
    structure = IStructure(lattice=lattice, species=["H", "He"], coords=coords)
    return SupercellInfo(structure,
                         "Pm-3m",
                         [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         sites)


def test_defect_entries_maker_from_sites(supercell_info):
    lattice = Lattice.cubic(1.0)
    coords = [[0.25, 0.25, 0.25]]
    structure = IStructure(lattice=lattice, species=["Li"], coords=coords)

    maker = ComplexDefectEntriesMaker(supercell_info,
                                      ComplexDefectSet([ComplexDefect("test", {0: None, 1: "Li"}, [0])]))
    assert maker.defect_entries == [DefectEntry(name="test",
                                                charge=0,
                                                structure=structure,
                                                site_symmetry="m-3m",
                                                defect_center=(0.125, 0.125, 0.125))]



"""
1. Make DefectEntriesMakeBySites
2. make a parser of defect_in_by_sites.yaml
3. make divacancy.in
"""