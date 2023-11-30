# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import re
from typing import Dict, Any

from pymatgen.core import Element


def remove_digits(name: str) -> str:
    """ "O1" -> "O" """
    return ''.join([i for i in name if not i.isdigit()])


def only_digits(name: str) -> str:
    """ "O1" -> "1" """
    return ''.join([i for i in name if i.isdigit()])


elements = [str(e) for e in Element]
e_va = elements + ["Va"]
e_i = elements + ["i"]


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


def typical_defect_name(name: str) -> bool:
    x = name.split("_")
    if len(x) == 2:
        _in, _out = x
        if _in in e_va and remove_digits(_out) in e_i:
            return True
    return False


def prettify_names(d: Dict[str, Any], style) -> Dict[str, Any]:
    result = {}
    out_names = [name.split("_")[1] for name in d.keys()]
    for name, v in d.items():
        in_name, out_name = name.split("_")
        r_out_name = remove_digits(out_name)
        out_name = r_out_name if f"{r_out_name}2" not in out_names else out_name
        _name = "_".join([in_name, out_name])
        if _name in result:
            raise ValueError("The prettified names are conflicted. "
                             "Change the defect names, please.")
        if style is None:
            pass
        elif style == "mpl":
            _name = defect_mpl_name(_name)
        else:
            raise ValueError(f"Style {style} is not adequate. Set mpl or None.")
        result[_name] = v
    return result


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


