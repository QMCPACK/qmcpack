.. _lab-qmc-basics:

Lab 2: QMC Basics
=================

Topics covered in this lab
--------------------------

This lab focuses on the basics of performing quality QMC calculations.
As an example, participants test an oxygen pseudopotential within DMC by
calculating atomic and dimer properties, a common step prior to
production runs. Topics covered include:

-  Converting pseudopotentials into QMCPACK’s FSATOM format

-  Generating orbitals with QE

-  Converting orbitals into QMCPACK’s ESHDF format with pw2qmcpack

-  Optimizing Jastrow factors with QMCPACK

-  Removing DMC time step errors via extrapolation

-  Automating QMC workflows with Nexus

-  Testing pseudopotentials for accuracy

Lab outline
-----------

#. Download and conversion of oxygen atom pseudopotential

#. DMC time step study of the neutral oxygen atom

   #. DFT orbital generation with QE

   #. Orbital conversion with

   #. Optimization of Jastrow correlation factor with QMCPACK

   #. DMC run with multiple time steps

#. DMC time step study of the first ionization potential of oxygen

   #. Repetition of a-d above for ionized oxygen atom

#. Automated DMC calculations of the oxygen dimer binding curve

Lab directories and files
-------------------------

::

  %
  labs/lab2_qmc_basics/
  │
  ├── oxygen_atom           - oxygen atom calculations
  │   ├── O.q0.dft.in          - Quantum ESPRESSO input for DFT run
  │   ├── O.q0.p2q.in          - pw2qmcpack.x input for orbital conversion run
  │   ├── O.q0.opt.in.xml      - QMCPACK input for Jastrow optimization run
  │   ├── O.q0.dmc.in.xml      - QMCPACK input file for neutral O DMC
  │   ├── ip_conv.py           - tool to fit oxygen IP vs timestep
  │   └── reference            - directory w/ completed runs
  │
  ├── oxygen_dimer          - oxygen dimer calculations
  │   ├── dimer_fit.py         - tool to fit dimer binding curve
  │   ├── O_dimer.py           - automation script for dimer calculations
  │   ├── pseudopotentials     - directory for pseudopotentials
  │   └── reference            - directory w/ completed runs
  │
  └── your_system           - performing calculations for an arbitrary system (yours)
      ├── example.py           - example nexus file for periodic diamond
      ├── pseudopotentials     - directory containing C pseudopotentials
      └── reference            - directory w/ completed runs

.. _lqb-pseudo:

Obtaining and converting a pseudopotential for oxygen
-----------------------------------------------------

First enter the ``oxygen_atom`` directory:

::

  cd labs/lab2_qmc_basics/oxygen_atom/

Throughout the rest of the lab, locations are specified with respect to ``labs/lab2_qmc_basics`` (e.g., ``oxygen_atom``).

We use a potential from the Burkatzki-Filippi-Dolg pseudopotential database.
Although the full database is available in QMCPACK distribution (``trunk/pseudopotentials/BFD/``),
we use a BFD pseudopotential to illustrate the process of converting and testing an
external potential for use with QMCPACK.   To obtain the pseudopotential, go to
http://www.burkatzki.com/pseudos/index.2.html
and click on the "Select Pseudopotential" button.  Next click on oxygen in the
periodic table.  Click on the empty circle next to "V5Z" (a large Gaussian
basis set) and click on "Next."  Select the Gamess format and click on
"Retrive Potential."  Helpful information about the pseudopotential will be
displayed.  The desired portion is at the bottom (the last 7 lines).  Copy
this text into the editor of your choice (e.g., ``emacs`` or ``vi``)
and save it as ``O.BFD.gamess``
(be sure to include a new line at the end of the file).  To transform the
pseudopotential into the FSATOM XML format used by QMCPACK, use the ``ppconvert``
tool:

::

  ppconvert --gamess_pot O.BFD.gamess --s_ref "1s(2)2p(4)" \
   --p_ref "1s(2)2p(4)" --d_ref "1s(2)2p(4)" --xml O.BFD.xml

Observe the notation used to describe the reference valence configuration for this helium-core PP: ``1s(2)2p(4)``.  The ``ppconvert`` tool uses the following convention for the valence states: the first $s$ state is labeled ``1s`` (``1s``, ``2s``, ``3s``, :math:`\ldots`), the first :math:`p` state is labeled ``2p`` (``2p``, ``3p``, :math:`\ldots`), and the first :math:`d` state is labeled ``3d`` (``3d``, ``4d``, :math:`\ldots`). Copy the resulting xml file into the ``oxygen_atom`` directory.

Note: The command to convert the PP into QE's UPF format is similar (both formats are required):

::

  ppconvert --gamess_pot O.BFD.gamess --s_ref "1s(2)2p(4)" \
   --p_ref "1s(2)2p(4)" --d_ref "1s(2)2p(4)" --log_grid --upf O.BFD.upf

For reference, the text of ``O.BFD.gamess`` should be:

::

  O-QMC GEN 2 1
  3
  6.00000000 1 9.29793903
  55.78763416 3 8.86492204
  -38.81978498 2 8.62925665
  1
  38.41914135 2 8.71924452

The full QMCPACK pseudopotential is also included in ``oxygen_atom/reference/O.BFD.*``.

.. _lqb-dft:

DFT with QE to obtain the orbital part of the wavefunction
----------------------------------------------------------

With the pseudopotential in hand, the next step toward a QMC calculation is to obtain the
Fermionic part of the wavefunction, in this case a single Slater determinant constructed
from DFT-LDA orbitals for a neutral oxygen atom.  If you had trouble with the pseudopotential conversion
step, preconverted pseudopotential files are located in the ``oxygen_atom/reference`` directory.

QE input for the DFT-LDA ground state of the neutral oxygen atom can be found in ``O.q0.dft.in``
and also in :ref:`Listing 58 <Listing 58>`.  Setting ``wf_collect=.true.`` instructs QE to write the
orbitals to disk at the end of the run. Option ``wf_collect=.true.`` could be a potential problem
in large simulations; therefore, we recommend avoiding it and using the converter pw2qmcpack in parallel
(see details in :ref:`pw2qmcpack`).
Note that the plane-wave energy cutoff has been set to a reasonable value of 300 Ry here (``ecutwfc=300``).
This value depends on the pseudopotentials used, and, in general,
should be selected by running DFT :math:`\rightarrow` (orbital conversion) :math:`\rightarrow` VMC with
increasing energy cutoffs until the lowest VMC total energy and variance is reached.

.. code-block::
  :caption: QE input file for the neutral oxygen atom (``O.q0.dft.in``)
  :name: Listing 58

  &CONTROL
     calculation       = 'scf'
     restart_mode      = 'from_scratch'
     prefix            = 'O.q0'
     outdir            = './'
     pseudo_dir        = './'
     disk_io           = 'low'
     wf_collect        = .true.
  /

  &SYSTEM
     celldm(1)         = 1.0
     ibrav             = 0
     nat               = 1
     ntyp              = 1
     nspin             = 2
     tot_charge        = 0
     tot_magnetization = 2
     input_dft         = 'lda'
     ecutwfc           = 300
     ecutrho           = 1200
     nosym             = .true.
     occupations       = 'smearing'
     smearing          = 'fermi-dirac'
     degauss           = 0.0001
  /

  &ELECTRONS
     diagonalization   = 'david'
     mixing_mode       = 'plain'
     mixing_beta       = 0.7
     conv_thr          = 1e-08
     electron_maxstep  = 1000
  /


  ATOMIC_SPECIES
     O  15.999 O.BFD.upf

  ATOMIC_POSITIONS alat
     O     9.44863067       9.44863161       9.44863255

  K_POINTS automatic
     1 1 1  0 0 0

  CELL_PARAMETERS cubic
          18.89726133       0.00000000       0.00000000
           0.00000000      18.89726133       0.00000000
           0.00000000       0.00000000      18.89726133

Run QE by typing

::

  mpirun -np 4 pw.x -input O.q0.dft.in >&O.q0.dft.out&

The DFT run should take a few minutes to complete.  If desired, you can track the progress of the DFT run by typing "``tail -f O.q0.dft.out``." Once finished, you should check the LDA total energy in ``O.q0.dft.out``  by typing  "``grep '!  ' O.q0.dft.out``."  The result should be close to

::

  !    total energy              =     -31.57553905 Ry

The orbitals have been written in a format native to QE in the ``O.q0.save`` directory.  We will convert them into the ESHDF format expected by QMCPACK by using the ``pw2qmcpack.x`` tool.  The input for ``pw2qmcpack.x`` can be found in the file ``O.q0.p2q.in`` and also in :ref:`Listing 59 <Listing 59>`.

.. code-block::
  :caption: ``pw2qmcpack.x`` input file for orbital conversion (``O.q0.p2q.in``)
  :name: Listing 59

  &inputpp
    prefix     = 'O.q0'
    outdir     = './'
    write_psir = .false.
  /

Perform the orbital conversion now by typing the following:

::

  mpirun -np 1 pw2qmcpack.x<O.q0.p2q.in>&O.q0.p2q.out&

Upon completion of the run, a new file should be present containing the orbitals for QMCPACK: ``O.q0.pwscf.h5``.  Template XML files for particle (``O.q0.ptcl.xml``) and wavefunction (``O.q0.wfs.xml``) inputs to QMCPACK should also be present.

.. _optimization-walkthrough:

Optimization with QMCPACK to obtain the correlated part of the wavefunction
---------------------------------------------------------------------------

The wavefunction we have obtained to this point corresponds to a noninteracting Hamiltonian.  Once the Coulomb pair potential is switched on between particles, it is known analytically that the exact wavefunction has cusps whenever two particles meet spatially and, in general, the electrons become correlated.  This is represented in the wavefunction by introducing a Jastrow factor containing at least pair correlations:

.. math::
  :label: eq66

  \Psi_{Slater-Jastrow}=e^{-J}\Psi_{Slater}

.. math::
  :label: eq67

  J = \sum_{\sigma\sigma'}\sum_{i<j}u^{\sigma\sigma'}_2(|r_i-r_j|) + \sum_\sigma\sum_{iI}u^{\sigma I}_1(|r_i-r_I|)\:.

Here :math:`\sigma` is a spin variable while :math:`r_i` and :math:`r_I`
represent electron and ion coordinates, respectively. The introduction
of :math:`J` into the wavefunction is similar to F12 methods in quantum
chemistry, though it has been present in essentially all QMC studies
since the first applications the method (circa 1965).

How are the functions :math:`u_2^{\sigma\sigma'}` and
:math:`u_1^{\sigma}` obtained? Generally, they are approximated by
analytical functions with several unknown parameters that are determined
by minimizing the energy or variance directly within VMC. This is
effective because the energy and variance reach a global minimum only
for the true ground state wavefunction
(:math:`\textrm{Energy}=E\equiv\langle{\Psi}|{\hat{H}}|{\Psi}\rangle`,
:math:`\textrm{Variance}=V\equiv\langle{\Psi}|{(\hat{H}-E)^2}|{\Psi}\rangle`).
For this exercise, we will focus on minimizing the variance.
