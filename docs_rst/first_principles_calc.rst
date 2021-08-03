Tips for first-principles calculations for point defects
--------------------------------------------------------

-----------------------------------------
1. How to treat symmetry of point defects
-----------------------------------------
As mentioned in the :doc:`tutorial`, the neighboring atoms near the defect are
initially slightly perturbed to break the symmetry.
However, some defects tend to move back to the symmetric atomic configuration
or recover a part of symmetry operations during the structure optimization.

Even in these cases, it is not apparent that the final structures are symmetric.
:code:`Pydefect` provides a sub-command that allows for symmetrizing the defect structure,
which is used as:

::

     pydefect_vasp_util rdp -p Va_Mg1_0/CONTCAR-finish -d Va_Mg1_0/defect_entry.json -n POSCAR_new

which create the POSCAR_new as follows.

::

    Mg7 Al16 O32
    1.0
    8.105607 0.000000 0.000000
    0.000000 8.105607 0.000000
    0.000000 0.000000 8.105607
    Mg Al O
    7 16 32
    direct
    0.749726 0.749738 0.249693 Mg
    0.749726 0.249738 0.749693 Mg
    0.249726 0.749738 0.749693 Mg
    0.003191 0.003203 0.003158 Mg
    0.496261 0.496273 0.003158 Mg
    0.496261 0.003203 0.496228 Mg
    0.003191 0.496273 0.496228 Mg
    0.625112 0.373467 0.373422 Al
    0.125997 0.874352 0.373422 Al
    0.125997 0.373467 0.874307 Al
    0.624659 0.874804 0.874760 Al
    0.624659 0.624671 0.624626 Al
    0.125997 0.126009 0.625079 Al
    0.125997 0.625124 0.125964 Al
    0.625112 0.126009 0.125964 Al
    0.373455 0.625124 0.373422 Al
    0.874340 0.126009 0.373422 Al
    0.874792 0.624671 0.874760 Al
    0.373455 0.126009 0.874307 Al
    0.373455 0.373467 0.625079 Al
    0.874792 0.874804 0.624626 Al
    0.874340 0.373467 0.125964 Al
    0.373455 0.874352 0.125964 Al
    0.862178 0.862190 0.862145 O
    0.361441 0.361453 0.859833 O
    0.361441 0.859878 0.361408 O
    0.859866 0.361453 0.361408 O
    0.388922 0.611788 0.611743 O
    0.887676 0.110542 0.611743 O
    0.887676 0.611788 0.110497 O
    0.393711 0.105752 0.105707 O
    0.611776 0.611788 0.388889 O
    0.105740 0.105752 0.393678 O
    0.110530 0.611788 0.887643 O
    0.611776 0.110542 0.887643 O
    0.611776 0.388934 0.611743 O
    0.110530 0.887688 0.611743 O
    0.105740 0.393723 0.105707 O
    0.611776 0.887688 0.110497 O
    0.637273 0.637285 0.862145 O
    0.138011 0.138023 0.859833 O
    0.138011 0.639597 0.361408 O
    0.639585 0.138023 0.361408 O
    0.637273 0.862190 0.637241 O
    0.138011 0.361453 0.639552 O
    0.138011 0.859878 0.137978 O
    0.639585 0.361453 0.137978 O
    0.393711 0.393723 0.393678 O
    0.887676 0.887688 0.388889 O
    0.887676 0.388934 0.887643 O
    0.388922 0.887688 0.887643 O
    0.862178 0.637285 0.637241 O
    0.361441 0.138023 0.639552 O
    0.361441 0.639597 0.137978 O
    0.859866 0.138023 0.137978 O

------------------------------------------
2. Tips for hybrid functional calculations
------------------------------------------
Hybrid functionals, especially the HSE06 functional,
and those with different exchange mixing parameters and/or screening distances,
have been regularly used for point-defect calculations recently.

Usually, hybrid functional calculations are a few tens more expensive
than those with a functional based on the local or semilocal density approximation.
Therefore, we need to take a little ingenuity to reduce their computational costs.

For this purpose, we regularly prepare the WAVECAR files obtained using a GGA functional.
(Although we can also relax the atomic positions using GGA beforehand,
it could be inappropriate for point-defect calculations,
because site symmetry of a defect calculated by GGA could be different from that by hybrid functionals.)

One can create the INCAR file for generating WAVECAR files using the GGA with the following command, for instance.

::

    grep -v LHFCALC INCAR | grep -v ALGO | sed s/"NSW     =  50"/"NSW     =   1"/ > INCAR-pre


