# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

import numpy as np
import pytest
from pydefect.analyzer.defect_structure_info import \
    Displacement, DefectStructureInfo, fold_coords, calc_drift, \
    calc_displacements, make_defect_structure_info, symmetrize_defect_structure
from pydefect.defaults import defaults
from pydefect.tests.helpers.assertion import assert_msonable
from pymatgen import Structure, Lattice
from vise.util.structure_symmetrizer import StructureSymmetrizer


@pytest.fixture
def structures():
    perfect = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Li", "U"],
        coords=[[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.5, 0.5, 0.5], [0, 0, 0]])

    initial = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.25, 0.25, 0.25], [0.75, 0.75, 0.75], [0.55, 0.5, 0.5], [0.0, 0.0, 0.0]])

    final = Structure(
        Lattice.cubic(10),
        species=["H", "He", "Be", "U"],
        coords=[[0.35, 0.35, 0.35], [0.74, 0.74, 0.74], [0.5, 0.5, 0.5], [0.0, 0.0, 0.0001]])
    return perfect, initial, final


@pytest.fixture
def displacements():
    x = np.sqrt((0.75-0.4)**2*3) * 10
    return [None,
            Displacement(specie="He",
                         original_pos=(0.75, 0.75, 0.75),
                         final_pos=(0.74, 0.74, 0.74),
                         distance_from_defect=round(x, 3),
                         displace_distance=round(np.sqrt(3) * 10 * 0.01, 3),
                         angle=0.0,
                         annotation="inward"),
            None,
            Displacement(specie="U",
                         original_pos=(0, 0, 0),
                         final_pos=(0.0, 0.0, 0.0001),
                         distance_from_defect=round(np.sqrt(3) * 10 * 0.4, 3),
                         displace_distance=0.001,
                         angle=54.7,
                         annotation=None),
            ]


@pytest.fixture
def def_str_info(structures, displacements):
    perfect, initial, final = structures
    return DefectStructureInfo(initial_point_group="m-3m",
                               final_point_group="m-3m",
                               initial_structure=initial,
                               final_structure=final,
                               perfect_structure=perfect,
#                               fd_to_id=[None, 1, 2, 3],
                               drift_dist=0.001,
#                               initial_defect_type=DefectType.substituted,
                               initial_vacancies=[2],
                               initial_interstitials=[2],
#                               final_defect_type=DefectType.complex,
                               final_vacancies=[0, 2],
                               final_interstitials=[0, 2],
                               defect_center_coord=(0.4, 0.4, 0.4),
                               displacements=displacements)
#                               neighboring_atom_indices=[0, 1, 2])


def test_fold_coords():
    actual = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                       species=["H"],
                       coords=[[0.5, 0, 0.49]])
    fold_coords(actual, [0, 0, 0])
    expected = Structure(Lattice.monoclinic(a=1, b=1, c=1, beta=20),
                         species=["H"],
                         coords=[[-0.5, 0, 0.49]])
    assert actual == expected


def test_calc_drift_dist(structures):
    perfect, _, final = structures
    drift_dist = calc_drift(perfect, final,
                            center=[0.4, 0.4, 0.4],
                            d_to_p=[None, 1, 2, 3])
    assert drift_dist == 0.001


def test_calc_displacements(structures, displacements):
    perfect, _, final = structures
    actual = calc_displacements(perfect, final,
                                center=[0.4, 0.4, 0.4],
                                d_to_p=[None, 1, None, 3])
    expected = displacements
    assert actual == expected


def test_str_info_msonable(def_str_info):
    assert_msonable(def_str_info)


def test_make_defect_structure_info(structures, def_str_info, mocker):
    perfect, initial, final = structures
    p_calc_results = mocker.stub(name='perfect_calc_results')
    defect_entry = mocker.stub(name='defect_entry')
    d_calc_results = mocker.stub(name='calc_results')

    p_calc_results.structure = perfect
    defect_entry.structure = initial
    defect_entry.site_symmetry = "m-3m"
    d_calc_results.structure = final
    d_calc_results.site_symmetry = "m-3m"
    actual = make_defect_structure_info(p_calc_results, defect_entry, d_calc_results)
    assert actual == def_str_info


def test_repr(def_str_info):
    expected = """
Site symmetry:
  m-3m -> m-3m
Defect center:
  [0.4, 0.4, 0.4]
Transferred atoms:
 - Cu [0.25, 0.25, 0.25]
 - Cu [0.25, 0.74, 0.74]
 + Cu [0.25, 0.5, 0.5]
 +  H [0.25, 0.25, 0.25]
Defect type:
  Complex defect
Sum of displacement distances:
  3.45 A
Displacements:
  Cu 1 [0.25, 0.75, 0.75] -> 1 [0.25, 0.75, 0.75] ([0.25, 0.75, 0.75])  0.0 A  0°
  Cu 2 [0.25, 0.75, 0.75] -> 1 [0.25, 0.75, 0.75] ([0.25, 0.75, 0.75])  0.0 A  0° inward
    """
    print(def_str_info)
#    assert def_str_info.__repr__() == expected


def test_symmetrize_defect_structure():
    structure = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
5 0 0
0 5 0
0 0 5
Mg   O
4     3
Direct
0.0051896248240553  0.9835077947659414  0.9945137498637422
0.0047282952713914  0.4827940046010823  0.4942929782542009
0.5040349492352973  0.9821499237428384  0.4944941755970405
0.5058945352747628  0.4828206016032297  0.9940420309511140
0.5045613848356609  0.4811103128264023  0.4933877756337353
0.0013796816599694  0.9829379087234287  0.4953360299212051
0.0083465288988691  0.4848714537370853  0.9941122597789658""")
    structure_symmetrizer = StructureSymmetrizer(
        structure,
        defaults.symmetry_length_tolerance,
        defaults.symmetry_angle_tolerance)
    actual = symmetrize_defect_structure(structure_symmetrizer,
                                         anchor_atom_idx=1,
                                         anchor_atom_coord=np.array([0.0, 0.5, 0.5]))
    expected = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
5 0 0
0 5 0
0 0 5
Mg   O
4     3
Direct
0.0 0.0 0.0
0.0 0.5 0.5
0.5 0.0 0.5
0.5 0.5 0.0
0.5 0.5 0.5
0.0 0.0 0.5
0.0 0.5 0.0""")
    assert actual == expected


def test_symmetrize_defect_structure_2():
    structure = Structure.from_str(fmt="POSCAR", input_string="""Mg48 N31
   1.00000000000000
    10.00    0.00    0.00
     0.00   10.00    0.00
     0.00    0.00   10.00
  N
    31
Direct
  0.9998735951099533  0.9999339271342436  0.9998810782738516
  0.4998325344890802  0.0000864604371742  0.0001093145892668
  0.9999579188430090  0.5001853666498661  0.0001277381372233
  0.2491639910909313  0.2182663238872422  0.4987861655656971
  0.2808293379596023  0.4991743721150215  0.7498359204125151
  0.5007987323114946  0.2498921962049820  0.2191539974347521
  0.2186260754052398  0.4998463318536253  0.2504951842089369
  0.5003683092799207  0.7505911171114192  0.2814698995089699
  0.7491029691281099  0.2824156531954642  0.4998653178588484
  0.2496769296641759  0.2810133141130393  0.0008972384265746
  0.2179934993920654  0.0013328906653882  0.7491830895564036
  0.9985995146190305  0.2494223137356002  0.2817274775328684
  0.2819549960242611  0.0002510594995684  0.2492863143591961
  0.9999066837513126  0.7494408251560003  0.2182162389686653
  0.7503446162775589  0.2186089947761758  0.0001821657373426
  0.5000178504783079  0.5000386610406409  0.9999895875233165
  0.4999380720704565  0.5002342427150381  0.5000689317878368
  0.0000976472392296  0.5000243131273407  0.5000777225283457
  0.5001616481443207  0.0002089601314523  0.4998675396277079
  0.7502599885437249  0.7191435719333086  0.9992528941462950
  0.7820064323950149  0.9990033992670391  0.2509026823008185
  0.0012722293791043  0.7506950871201497  0.7182763220765622
  0.7176368346430237  0.9998582962107960  0.7509680009789932
  0.0000228430868177  0.2509182464808006  0.7821761073165092
  0.2495215811710665  0.7814963684974714  0.9998566240987685
  0.7508300518084354  0.7818602560717594  0.5013867902350881
  0.7190618878688895  0.5010405127949369  0.2502514755283229
  0.4989978969018409  0.7502977850544852  0.7809492219327865
  0.7814464623477875  0.5003886730650109  0.7494947707104060
  0.4996606931909255  0.2496508616713697  0.7186036919929819
  0.2506716727065808  0.7181216545822622  0.5001902272634595""")
    structure_symmetrizer = StructureSymmetrizer(
        structure,
        defaults.symmetry_length_tolerance,
        defaults.symmetry_angle_tolerance)
    actual = symmetrize_defect_structure(structure_symmetrizer,
                                         anchor_atom_idx=15,
                                         anchor_atom_coord=np.array([0.5, 0.5, 0.0]))

    expected = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
10 0 0
0 10 0
0 0 10
N
31
Direct
    0         0         0
    0.5       0         0
    0         0.5       0
    0.249164  0.218235  0.498835
    0.280829  0.499143  0.749885
    0.500857  0.249885  0.219171
    0.218626  0.499815  0.250544
    0.500185  0.750544  0.281374
    0.749103  0.282384  0.499914
    0.249885  0.280829  0.000857
    0.218235  0.001165  0.749164
    0.998835  0.249164  0.281765
    0.282384  8.6e-05   0.249103
    0.999914  0.749103  0.217616
    0.750544  0.218626  0.000185
    0.5       0.5       0
    0.5       0.5       0.5
    0         0.5       0.5
    0.5       1         0.5
    0.750115  0.719171  0.999143
    0.781765  0.998835  0.250836
    0.001165  0.750836  0.718235
    0.717616  0.999914  0.750897
    8.6e-05   0.250897  0.782384
    0.249456  0.781374  0.999815
    0.750836  0.781765  0.501165
    0.719171  0.500857  0.250115
    0.499143  0.750115  0.780829
    0.781374  0.500185  0.749456
    0.499815  0.249456  0.718626
    0.250897  0.717616  0.500086""")
    assert actual == expected


def test_symmetrize_defect_structure_wo_anchor():
    structure = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
5 0 0
0 5 0
0 0 5
Mg   O
4     3
Direct
0.01 0.01 0.01
0.01 0.51 0.51
0.51 0.01 0.51
0.51 0.51 0.01
0.51 0.51 0.51
0.01 0.01 0.51
0.01 0.51 0.01""")
    structure_symmetrizer = StructureSymmetrizer(structure)
    actual = symmetrize_defect_structure(structure_symmetrizer=structure_symmetrizer)
    expected = Structure.from_str(fmt="POSCAR", input_string="""Mg4 O3
1.00000000000000
5 0 0
0 5 0
0 0 5
Mg   O
4     3
Direct
0.01 0.01 0.01
0.01 0.51 0.51
0.51 0.01 0.51
0.51 0.51 0.01
0.51 0.51 0.51
0.01 0.01 0.51
0.01 0.51 0.01""")
    assert actual == expected