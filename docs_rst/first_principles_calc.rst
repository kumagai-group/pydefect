Tips for first-principles calculations for point defects
--------------------------------------------------------

As mentioned in the :doc:`tutorial`, the neighboring atoms near the defect
are initially slightly perturbed to break the symmetry.
And, some defects tend to move back to the symmetric atomic configuration during the structure optimization.
Some other defects might recover a part of symmetry operations.

Even in these cases, the s



Tips for hybrid functional calculations
---------------------------------------

Hybrid functionals, especially HSE06 functional
and those with different exchange mixing parameters and/or screening distances,
have been regularly used for point-defect calculations.

Usually, calculations with a hybrid functional are a few tens more expensive
than those with a functional based on the local or semilocal density approximation.
Therefore, we need to take a little ingenuity to reduce their computational costs.

For this purpose, we regularly prepare the WAVECAR files obtained using GGA.
(Although we also relax the atomic positions using GGA beforehand in some cases, it would be inappropriate for point-defect
calculations, because site symmetry of a defect calculated by GGA could be different from that by hybrid functionals.
Furthermore, electronic structures of defects could also be different.)

One can create the INCAR file for generating WAVECAR files using the GGA with the following command, for instance.

::
    grep -v LHFCALC INCAR | grep -v ALGO | sed s/"NSW     =  50"/"NSW     =   1"/ > INCAR-pre


