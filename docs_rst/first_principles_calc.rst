Tips for first-principles calculations for point defects
--------------------------------------------------------

----------------------------------
1. Symmetrization of point defects
----------------------------------
As mentioned in the :doc:`tutorial`, the atoms neighboring a defect
are initially slightly displaced to break the symmetry.
However, some defects tend to relax back to a symmetric configuration,
or partially recover symmetry operations during structure optimization.

Even in such cases, the relaxed structures are not obviously symmetric.
:code:`Pydefect` provides the :code:`refine_defect_poscar` (= :code:`rdp`) sub-command,
which enables symmetrization of the defect structure.
It is used as follows:

::

     pydefect_vasp_util rdp -p Va_Mg1_0/CONTCAR-finish -d Va_Mg1_0/defect_entry.json -n POSCAR_new

This command creates the :code:`POSCAR_new` file, as shown below:

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
    ...
    (truncated for brevity)
    ...

------------------------------------------
2. Tips for hybrid functional calculations
------------------------------------------
Hybrid functionals—such as the HSE06 functional—
and their variants with different exchange mixing parameters or screening lengths
have recently become standard tools for point-defect calculations.

These calculations are typically tens of times more computationally expensive
than those based on local or semilocal density functional approximations.
Thus, some ingenuity is needed to reduce their computational cost.

To address this, it is common to prepare a :code:`WAVECAR` file
from a GGA-based calculation beforehand.
(Note that relaxing atomic positions using GGA may be inappropriate for point-defect calculations,
as the site symmetry of a defect calculated by GGA can differ from that of a hybrid functional.
If the symmetry is increased by GGA, a hybrid functional calculation will not reduce it again.)

You can generate an :code:`INCAR` file for creating the :code:`WAVECAR`
using GGA with the following command:

::

    grep -v LHFCALC INCAR | grep -v ALGO | sed s/"NSW     =  50"/"NSW     =   1"/ > INCAR-pre

Then, move the resulting :code:`WAVECAR` file to the directory
where the hybrid functional calculation will be performed.
