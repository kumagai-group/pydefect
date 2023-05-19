# # -*- coding: utf-8 -*-
# #  Copyright (c) 2022 Kumagai group.
# import pytest
# from pydefect.analyzer.defect_energy import CrossPoints
# from pydefect.analyzer.transition_levels import make_transition_levels, \
#     show_transition_levels, TransitionLevel, TransitionLevels
# from pydefect.analyzer.transition_levels_plotter import \
#     TransitionLevelsMplPlotter, make_mpl_tl_data, MplTLData


# def test_make_mpl_tl_data():
#     tl1 = TransitionLevel("Va_O1", [[2, 1], [1, 0]], [2.0, 3.0], [1.0, 2.0])
#     tl2 = TransitionLevel("Va_Mg1", [[-2, -1]], [2.0], [1.0])
#     tls = TransitionLevels([tl1, tl2], cbm=3.0, supercell_vbm=-0.5, supercell_cbm=0.5)
#     actual = make_mpl_tl_data(tls)
#     expected = MplTLData(labels=["Va_O1", "Va_Mg1"],
#                          charges=[2, 1, 0, -1, -2],
#                          tl_energy_widths=[[1.0, 0.0], [1.0, 0.0], [1.0, 0.0],
#                                            [0.0, 2.0], [0.0, 1.0]],
#                          piled_tl_energy_widths=[[1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [3.0, 2.0], [3.0, 3.0]])
#     assert actual == expected

# # def test_plot_transition_levels():
# #     # tl1 = TransitionLevel("Va_O1", [[2, 1], [1, 0]],
# #     #                       [1.23456789, 2.23456789], [3.23456789, 4.23456789])
# #     # tl2 = TransitionLevel("Va_Mg1", [[-2, -1], [-1, 0]],
# #     #                       [10.2345678, 20.2345678], [30.2345678, 40.2345678])

#     # # actual = plot_transition_levels([tl1, tl2])
#     # TransitionLevelsMplPlotter().plt.show()