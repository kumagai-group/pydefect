Tutorial of how to calculate carrier and defect concentration using pydefect
--------------------

Theory how to calculate the carrier and defect concentrations
is written in our paper `[Y. Kumagai et al., Phys. Rev. B 90, 125202 (2014)]  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.125202>`_.


=========================
1. Carrier concentrations
=========================
In semiconductors, the Fermi level, :math:`E_mathrm{F}`, normally lies inside the
band gap.  The electron (:math:`n`) and hole (:math:`p`) concentrations are
related to :math:`E_\mathrm{F}` through the Fermi–Dirac distribution:

.. math::

    n = \int_{\epsilon_\text{CBM}}^{\infty}
        D(\epsilon)\,
        \frac{1}{e^{(\epsilon - E_\mathrm{F})/k_B T} + 1}\, d\epsilon ,

and

.. math::

    p = \int_{-\infty}^{\epsilon_\text{VBM}}
        D(\epsilon)\,
        \frac{1}{e^{(E_\mathrm{F} - \epsilon)/k_B T} + 1}\, d\epsilon ,

where :math:`D(\epsilon)` is the electronic density of states (DOS) and
:math:`T` is the temperature.  At a fixed :math:`T`, the Fermi level therefore
uniquely determines :math:`n` and :math:`p`.

For their calculations, the material-dependent quantity is only the DOS, and
thus, initially, we need to calculate it.

After performing the DOS calculations,
using the following command,

::

    pydefect_vasp_util make_total_dos ../unitcell/dos/vasprun.xml -o ../unitcell/dos/OUTCAR

we can create a :code:`total_dos.json` file.
Then, using the following command,

::

    pydefect_util calc_carrier_concentrations -t total_dos.json  -T 1000

we can calculate the carrier concentrations,
which is written in :code:`con_by_Ef_only_pn_1000.0K.json`.
Usually, we plot at the synthesized temperature and room temperature, namely 300K.
These data can be plotted using

::

    pydefect_util plot_carrier_concentrations

========================
2. Defect concentrations
========================
We then explain how to calculate the defect concentrations.

Since the defect concentrations are related to the
number of equivalent sites, which are also linked with the site symmetry,
initially, we evaluate the site degeneracies in the supercell using

::

    pydefect_util make_degeneracies -d ../defects/*_*_*/ -s ../defects/supercell_info.json

With this command, we can create :code:`degeneracies.yaml`.

Next, we can calculate the defect concentration at a position in the chemical potential diagram
using this command

::

    pydefect_util calc_defect_concentrations -d ../defects/defect_energy_summary.json -l A --degeneracies degeneracies.yaml -t total_dos.json

This calculation requires iterative calculation as follows.

::

    INFO: Calc equilibrium concentration: iteration 0
    INFO: Calc equilibrium concentration: iteration 1
    INFO: Calc equilibrium concentration: iteration 2
    INFO: Calc equilibrium concentration: iteration 3
    INFO: Calc equilibrium concentration: iteration 4
    INFO: Calc equilibrium concentration: iteration 5
    INFO: Calc equilibrium concentration: iteration 6
    INFO: Equilibrium concentration is found.
    T: 300
    Pinning levels:
    -----------  -----  -----
    Cs_2VCu             3.150
    Cs_Cu1
    Cs_i1        1.266
    Cu_Cs1
    Cu_i1
    O_i1
    Va_Cs1              3.201
    Va_Cu-split
    Va_Cu1              2.796
    Va_O1        0.028
    -----------  -----  -----

In addition, we can calculate the carrier concentrations at 300~K
on the premise that defect concentrations are quenched from 1000~K,
while the defects can occupy different charge states in the ratio determined by the Boltzmann distribution.
For this, we can use :code:`--con_by_Ef` option.

==========================================
3. Technical advice about DFT calculations
==========================================

In defect calculations, we usually confine :math:`E_\mathrm{F}` to the
calculated band gap.  Two situations, however, require us to place
:math:`E_\mathrm{F}` outside the gap.

* **Degenerate doping.**
  When charged defects drive :math:`E_\mathrm{F}` into a band, the
  semiconductor becomes *degenerate*.  For nitride semiconductors such as InN,
  Zn\ :math:`_3`\ N\ :math:`_2`, and ScN the positively charged nitrogen
  vacancy forms easily, donates many electrons, and pushes
  :math:`E_\mathrm{F}` well above the conduction-band minimum (CBM).  The
  resulting band filling produces the Burstein–Moss (BM) blue-shift, so that
  optical gaps exceed the fundamental gap.

* **Band-edge correction.**
  Semi-local DFT functionals (e.g. PBE) underestimate the band gap.
  A common remedy is to rigidly shift the valence- and conduction-band edges to
  match more accurate hybrid-DFT or :math:`G_0W_0` gaps.  After the shift
  :math:`E_\mathrm{F}` may lie above the corrected CBM or below the corrected
  valence-band maximum (VBM).

When computing defect formation energies and thermodynamic transition levels
(TTLs) for such cases, the :math:`k`-point mesh needs careful attention.  If
the mesh contains the true band extremum (e.g. the :math:`\Gamma` point for a
:math:`\Gamma`-centred CBM), an electronic level that should sit above the CBM
can relax into the CBM of the supercell.  The released electrons occupy a
perturbed host state (PHS), pinning the TTL artificially close to the CBM.

A practical workaround is to use a Monkhorst–Pack grid that **excludes** the
band extremum.  For instance, a :math:`2 \times 2 \times 2` mesh omitting the
:math:`\Gamma` point widens the supercell gap spuriously, but also keeps donor
TTLs inside this artificial gap.  Defect formation energies evaluated with
such a mesh remain physically meaningful when
:math:`E_\mathrm{F}` lies outside the true band edges.
