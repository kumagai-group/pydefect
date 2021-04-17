# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
from pymatgen import Element


def remove_digits(name: str) -> str:
    """ "O1" -> "O" """
    return ''.join([i for i in name if not i.isdigit()])


def only_digits(name: str) -> str:
    """ "O1" -> "1" """
    return ''.join([i for i in name if i.isdigit()])


elements = [str(e) for e in Element]


def defect_mpl_name(name: str) -> str:
    """ "Va_O1" -> "$V_{{\\rm O}1}$"
        "Mg_i1" -> "${\\rm Mg}_{i1}$" """
    in_name, out_name = name.split("_")
    if in_name in elements:
        in_name = "{\\rm " + in_name + "}"
    elif in_name == "Va":
        in_name = "V"

    r_out_name = remove_digits(out_name)
    if r_out_name in elements:
        out_name = "{{\\rm " + r_out_name + "}" + only_digits(out_name) + "}"
    else:
        out_name = "{" + out_name + "}"

    return f"${in_name}_{out_name}$"


def defect_plotly_name(name: str) -> str:
    """ "Va_O1" -> "<i>V</i><sub>O1</sub>" """
    in_name, out_name = name.split("_")
    if in_name == "Va":
        in_name = "<i>V</i>"
    out_name = f"<sub>{out_name}</sub>"
    return f"{in_name}{out_name}"


def defect_plotly_full_name(fullname: str) -> str:
    """ "Va_O1_1" -> "<i>V</i><sub>O1</sub><sup>1</sup>" """
    name, charge = fullname.rsplit('_', 1)
    return f"{defect_plotly_name(name)}<sup>{charge}</sup>"


# def defect_html_title_name(fullname):
#     x = fullname.split("_")
#     if len(x) == 2:
#         in_name, out_name = x
#     elif len(x) == 3:
#         in_name, out_name, charge = x
#     else:
#         raise ValueError
#
#     if in_name == "Va":
#         in_name = html.I("V")
#     else:
#         in_name = html.Span(in_name)
#
#     result = [in_name, html.Sub(out_name)]
#     if len(x) == 3:
#         result.append(html.Sup(charge))
#     return result


# def sanitize_defect_energies_for_plot(defect_energies: List[DefectEnergy],
#                                       for_plotly: bool = False):
#     result = []
#     out_names = [e.name.split("_")[1] for e in defect_energies]

# for e in defect_energies:
#     ee = deepcopy(e)
#     in_name, out_name = e.name.split("_")
#     r_out_name = remove_digits(out_name)
#     out_name = r_out_name if f"{r_out_name}2" not in out_names else out_name
#     if for_plotly:
#         ee.name = defect_plotly_name("_".join([in_name, out_name]))
#     else:
#         ee.name = defect_mpl_name("_".join([in_name, out_name]))
#     result.append(ee)

# return result

