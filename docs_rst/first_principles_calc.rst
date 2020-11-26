Tips for first-principles calculations for point defects
--------------------------------------------------------

-----------------------------------------
1. How to treat symmetry of point defects
-----------------------------------------
As mentioned in the :doc:`tutorial`, the neighboring atoms near the defect are initially slightly perturbed to break the symmetry.
And, some defects tend to move back to the symmetric atomic configuration during the structure optimization.
Some other defects might recover a part of symmetry operations.

Even in these cases, the final structures do not look symmetric.
When one wants to symmetrize the defect structure, use the following command.

::

    python $PYDEFECT_PATH/pydefect/cli/vasp/util_commands/make_refined_poscar.py

which will create the symmetrized :code:`POSCAR` file if the structure is not the P1 symmetry.
Then, the previous :code:`OUTCAR` and :code:`CONTCAR` are moved to
:code:`OUTCAR.sym_1` and:code:`CONTCAR.sym_1`, respectively.

It is also possible to include the command in the runshell script as an example:

::

    $VASP_cmd

    hostname > host
    name=`basename "$PWD"`
    echo $name
    if [ $name != "perfect" ]; then
        python $PYDEFECT_PATH/pydefect/cli/vasp/util_commands/make_refined_poscar.py
        if [ -e CONTCAR.sym_1 ]; then
            $VASP_cmd
        fi
    fi

------------------------------------------
2. Tips for hybrid functional calculations
------------------------------------------
Hybrid functionals, especially the HSE06 functional
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


