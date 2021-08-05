Tutorial for calculation of vertical transition level
-----------------------------------------------------

We here illustrate how to calculate the vertical transition level (VTL)
with an example of an oxygen vacancy in MgAl2O4.
For the calculation of the VTL, we need to apply special correction scheme,
here we call it the GKFO correction.
Please read
`T. Gake, Y. Kumagai*, C. Freysoldt, and F. Oba, Phys. Rev. B, 101, 020102(R) (2020).
<link.aps.org/doi/10.1103/PhysRevB.101.020102>`_
for details.

Assuming that the defect calculations in MgAl2O4 based on the PBEsol functional
have been done already as introduced in the tutorial
and one further wants to calculate the photo-absorption energy caused by excitation
of an oxygen vacancy localized state in the neutral charge state to the CBM.

We here assume the following directory structure.

::

    MgAl2O4
     │
     ├ unitcell/ ── unitcell.json
     │
     └ defects/ ── perfect/
                 └ Va_O1_0/ ── absorption/

Firstly, we create the :code:`absorption/` directory at the :code:`Va_O1_0/`
and copy the vasp input files from :code:`Va_O1_0/`.
Then, edit :code:`INCAR` to change :code:`NSW` to 1 and add the :code:`NELECT` tag with
reducing it by 1 from the neutral one.
Or one can type:

::

    vise vs -uis NSW 1 --options charge 1 -d ../ -t defect

After running vasp, we then create :code:`calc_results.json`
using the following command in the :code:`absorption/` directory.

::

    pydefect cr -d .


And, wee obtain :code:`gkfo_correction.pdf` and :code:`gkfo_correction.json` files with the following command.

::

    pydefect_util gkfo -u ../../../unitcell/unitcell.yaml -iefnv ../correction.json -icr ../calc_results.json -fcr calc_results.json -cd 1

With the :code:`pydefect_print` command, the correction energy is shown as follows:


::

    +--------------------+------------+
    | charge             |  0         |
    | additional charge  |  1         |
    | pc 1st term        |  0         |
    | pc 2nd term        |  0.819326  |
    | alignment 1st term | -0.0953654 |
    | alignment 2nd term | -0.10711   |
    | alignment 3rd term | -0         |
    | correction energy  |  0.61685   |
    +--------------------+------------+


The :code:`gkfo_correction.pdf` shows the potential profile caused by addition/removal
of the electron and its alignment term.

.. image:: gkfo_correction.pdf


For the absorption energy, one needs to know the conduction band minimum position,
which is now 9.2376 eV in MgAl2O4. And, the total energies of initial and final states are
-405.739 and -411.920 eV, respectively.
Therefore, the absorption energy is

::

    -411.920 + 405.739 + 9.2376 + 0.61685 =  3.67345 eV

It is also worthwhile to check the eigenvalues of initial and final states.

We can obtain the :code:`eigenvalues.pdf`, which looks as

.. image:: eigenvalues_final.pdf

and the initial :code:`eigenvalues.pdf` looks as

.. image:: eigenvalues_initial.pdf


