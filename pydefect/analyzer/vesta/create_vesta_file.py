# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
# This file is originally developed by Naoki Tsunoda.

import itertools
from typing import Union, Optional, Dict, Iterable

from pymatgen.core.structure import Structure
from vise.util.logger import get_logger

logger = get_logger(__name__)


def replace_dummy_to_xx(target_str: str) -> str:
    return target_str.replace("X0+", "XX", -1)


def val_to_str_line(val_line: list or tuple) -> str:
    format_str = "{{:.{0}f}}".format(6)
    return " ".join([format_str.format(c) for c in val_line])


class VestaFile:
    def __init__(self,
                 structure: Structure,
                 title: str = None,
                 vectors: Dict[int, Iterable] = None,
                 boundary: Optional[Iterable] = None):
        """
        Args:
            vectors (dict):
                e.g., {1: [1.0, 0.0, 0.0], ..}
            boundary (Iterable):
                a_min, a_max, b_min, b_max, c_min, c_max
                Center of origin is determined from the mean of boundary.
                e.g., "-0.5 0.5 -0.5 0.5 -0.5 0.5" means (0, 0, 0) as center.
        """
        self.blocks = [Title(structure, title),
                       Cellp(structure),
                       Struc(structure),
                       Bound(boundary),
                       SBond(structure),
#                       DummyAtomt(structure),
                       Vect(vectors),
                       Style(structure)]

    def __repr__(self):
        outs = []
        for block in self.blocks:
            outs.append(repr(block))
        return "\n".join(outs)

    def write_file(self, filename: str):
        filename = filename if ".vesta" in filename else f"{filename}.vesta"
        with open(filename, 'w') as poscar_vesta:
            poscar_vesta.write(repr(self))
        logger.info(f"{filename} is created.")


class Title:
    """ TITLE block in *.vesta files showing title in VESTA. """
    vesta_version = "#VESTA_FORMAT_VERSION 3.5.0"
    header = "TITLE"

    def __init__(self, structure: Structure, title: str = None):
        self.title = title or structure.formula

    def __repr__(self):
        return "\n".join([f'{self.vesta_version}',
                          f'{self.header}',
                          f'{self.title}'])


class Cellp:
    """
    CELLP block in *.vesta files, describing lattice constants.
    # e.g.) 3.992878   3.992878   3.992878  90.000000  90.000000  90.000000
    """
    header = "CELLP"

    def __init__(self, structure: Structure):
        self.cell_param = structure.lattice.lengths + structure.lattice.angles

    def __repr__(self):
        return "\n".join([self.header, val_to_str_line(self.cell_param)])


class Struc:
    """
    STRUC block in *.vesta files, showing atom sites.
    The index starts from 1
    e.g. "index" "element" "element+index" + "occupation" + "coords"
          1       Ba        Ba1               1.0000         0.500  0.500  0.500
          0 0 0 0 0
    """
    header = "STRUC"
    separator = " 0 0 0 0 0 "
    zero_coord = " 0.0 0.0 0.0 "

    def __init__(self, structure: Structure):
        coords = []
        for site, c in enumerate(structure, 1):
            specifies = f'{site} {c.species_string} ' \
                        f'{c.species_string}{site} {1.0} '
            coord = val_to_str_line(c.frac_coords)
            coords.append(specifies + coord)
            coords.append(self.zero_coord)
        # replace "X0+" (pmg dummy species) -> "XX" (vesta dummy species)
        self.str_coords = replace_dummy_to_xx('\n'.join(coords))

    def __repr__(self):
        outs = [self.header, self.str_coords, self.separator]
        return '\n'.join(outs)


class Bound:
    """control range of plot"""
    header = "BOUND"
    separator = " 0 0 0 0 0 "

    def __init__(self,
                 boundary: Optional[Union[tuple, list]] = (0, 1, 0, 1, 0, 1)):
        """
        Args:
            boundary (tuple or list): a_min, a_max, b_min, b_max, c_min, c_max
        """
        if boundary:
            if len(boundary) != 6:
                raise ValueError("length of boundary must be 6")
            self.boundary = val_to_str_line(boundary)
        else:
            self.boundary = None

    def __repr__(self):
        if self.boundary:
            outs = [self.header, self.boundary, self.separator]
            return "\n".join(outs)
        else:
            return ""


class SBond:
    """
    SBOND block in *.vesta files, control show/hide bonding between elements.
    """
    header = "SBOND"
    separator = " 0 0 0 0 "

    def __init__(self, structure: Structure, bond_factor=1.2):
        bond_list = []
        i = 1
        for e1, e2 in itertools.permutations(structure.types_of_species, 2):
            try:
                bond = float(e1.average_ionic_radius + e2.average_ionic_radius) * bond_factor
            except AttributeError:
                continue
            bond_list.append(
                f"{i} {e1.symbol} {e2.symbol} 0.0 {bond:5.4}  0  1  1  0  1")
            i += 1
        self.bonds = "\n".join(bond_list)

    def __repr__(self):
        outs = [self.header, self.bonds, self.separator]
        return "\n".join(outs)


class Vect:
    """
     VECTR block in *.vesta files, showing directions of vectors
    e.g)
    1(vct_index)   -0.16519    0.00000    0.00000  <-- arrow end point
    1(atm_index)   0           0          0        <-- options

            penetration add_atomic_radius
    type 0:    no            no
    type 1:   yes            no
    type 2:    no           yes
    type 3:   yes           yes

    """
    header_1 = "VECTR"
    separator = " 0 0 0 0 0 "
    header_2 = "VECTT"
    type = 2

    def __init__(self,
                 vectors: Dict[int, Iterable],
                 size: float = 0.5,
                 color: str = "1 1 1"):
        """
        Args:
            vectors (dict): dict of site index and vector.
                            e.g., {1: [1.0, 0.0, 0.0], ..}
        """
        if vectors:
            vec_lines = []
            vector_types = []
            vector_option = f"{size} {color} {self.type}"
            for i, (idx, vec) in enumerate(vectors.items(), 1):
                vec_lines.append(f'{i} {val_to_str_line(vec)}')
                vec_lines.append(f"{idx}  0 0 0 0")
                vec_lines.append(self.separator)
                vector_types.append(f'{i} {vector_option}')
            self.vectors = '\n'.join(vec_lines)
            self.vector_types = '\n'.join(vector_types)
        else:
            self.vectors = None

    def __repr__(self):
        if self.vectors:
            return "\n".join([self.header_1, self.vectors, "",
                              self.header_2, self.vector_types])
        else:
            return ""


class DummyAtomt:
    """
    This is class object of ATOMT block in *.vesta files.
    Control style of atom object.
    (e.g., 1         Ba  0.8000  30 239  44  30 239  44 204)
    """
    header = "ATOMT"

    def __init__(self, structure: Structure):
        if "X" in structure.symbol_set:
            self.str = "  1         Xx  0.2000  76  76  76  76  76  76 504"
        else:
            self.str = ""

    def __repr__(self):
        return "\n".join([self.header, self.str]) if self.str else ""


class Style:
    """
    STYLE block in *.vesta files.
    """
    header = "STYLE"
    amp_prefix = "VECTS "
    atoms_prefix = "ATOMS "
    sects_prefix = "SECTS"
    sectp_prefix = "SECTP"
    ucolp_prefix = "UCOLP"

    all_bold_cell = "0  2  1.000   0   0   0"
    isurf = "1   1    12.0991 255 255   0 127 255"
    sects = " 160  1"
    sectp = "  1 0 0  0.00000E+00  0.00000E+00  0.00000E+00  0.00100E+00"

    def __init__(self, structure: Structure, is_ionic: bool = True):
        self.amplitude = (structure.volume / structure.num_sites) ** 1 / 3
        self.atoms = "1  0  1" if is_ionic else "0  0  1"

    def __repr__(self):
        vct = f'{self.amp_prefix} {self.amplitude}'
        atom = f'{self.atoms_prefix} {self.atoms}'
        sects = f'{self.sects_prefix} {self.sects}'
        sectp = f'{self.sectp_prefix} \n {self.sectp}'
        ucol = f'{self.ucolp_prefix} \n {self.all_bold_cell}'
        return "\n".join([self.header, vct, sects, sectp, ucol, atom])
