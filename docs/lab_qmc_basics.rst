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

First, we need to update the template particle and wavefunction information in ``O.q0.ptcl.xml`` and ``O.q0.wfs.xml``.  We want to simulate the O atom in open boundary conditions (the default is periodic).  To do this, open ```O.q0.ptcl.xml`` with your favorite text editor (e.g., ``emacs`` or ``vi``) and replace

::

  <parameter name="bconds">
     p p p
  </parameter>
  <parameter name="LR_dim_cutoff">
     15
  </parameter>

with

::

  <parameter name="bconds">
     n n n
  </parameter>

Next we will select Jastrow factors appropriate for an atom.  In open boundary conditions, the B-spline Jastrow correlation functions should cut off to zero at some distance away from the atom.  Open ``O.q0.wfs.xml`` and add the following cutoffs (``rcut`` in Bohr radii) to the correlation factors:

::

  ...
  <correlation speciesA="u" speciesB="u" size="8" rcut="10.0">
  ...
  <correlation speciesA="u" speciesB="d" size="8" rcut="10.0">
  ...
  <correlation elementType="O" size="8" rcut="5.0">
  ...

These terms correspond to
:math:`u_2^{\uparrow\uparrow}/u_2^{\downarrow\downarrow}`,
:math:`u_2^{\uparrow\downarrow}`, and
:math:`u_1^{\uparrow O}/u_1^{\downarrow O}`, respectively. In each case,
the correlation function (:math:`u_*`) is represented by piecewise
continuous cubic B-splines. Each correlation function has eight
parameters, which are just the values of :math:`u` on a uniformly spaced
grid up to ``rcut``. Initially the parameters (``coefficients``) are set
to zero:

::

  <correlation speciesA="u" speciesB="u" size="8" rcut="10.0">
    <coefficients id="uu" type="Array">
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    </coefficients>
  </correlation>

Finally, we need to assemble particle, wavefunction, and pseudopotential information into the main QMCPACK input file (``O.q0.opt.in.xml``) and specify inputs for the Jastrow optimization process.  Open ``O.q0.opt.in.xml`` and write in the location of the particle, wavefunction, and pseudopotential files ("``<!-- ... -->``" are comments):

::

  ...
  <!-- include simulationcell and particle information from pw2qmcpqack -->
  <include href="O.q0.ptcl.xml"/>
  ...
  <!-- include wavefunction information from pw2qmcpqack -->
  <include href="O.q0.wfs.xml"/>
  ...
  <!-- O pseudopotential read from "O.BFD.xml" -->
  <pseudo elementType="O" href="O.BFD.xml"/>
  ...

The relevant portion of the input describing the linear optimization process is

::

  <loop max="MAX">
    <qmc method="linear" move="pbyp" checkpoint="-1">
      <cost name="energy"              >  ECOST    </cost>
      <cost name="unreweightedvariance">  UVCOST   </cost>
      <cost name="reweightedvariance"  >  RVCOST   </cost>
      <parameter name="timestep"       >  TS       </parameter>
      <parameter name="samples"        >  SAMPLES  </parameter>
      <parameter name="warmupSteps"    >  50       </parameter>
      <parameter name="blocks"         >  200      </parameter>
      <parameter name="subSteps"       >  1        </parameter>
      <parameter name="nonlocalpp"     >  yes      </parameter>
      <parameter name="useBuffer"      >  yes      </parameter>
      ...
    </qmc>
  </loop>

An explanation of each input variable follows.  The remaining variables control specialized internal details of the linear optimization algorithm.  The meaning of these inputs is beyond the scope of this lab, and reasonable results are often obtained keeping these values fixed.

energy
   Fraction of trial energy in the cost function.

unreweightedvariance
   Fraction of unreweighted trial variance in the cost function.
   Neglecting the weights can be more robust.

reweightedvariance
   Fraction of trial variance (including the full weights) in the cost
   function.

timestep
   Time step of the VMC random walk, determines spatial distance moved
   by each electron during MC steps. Should be chosen such that the
   acceptance ratio of MC moves is around 50% (30–70% is often
   acceptable). Reasonable values are often between 0.2 and 0.6
   :math:`\textrm{Ha}^{-1}`.

samples
   Total number of MC samples collected for optimization; determines
   statistical error bar of cost function. It is often efficient to
   start with a modest number of samples (50k) and then increase as
   needed. More samples may be required if the wavefunction contains a
   large number of variational parameters. MUST be be a multiple of the
   number of threads/cores.

warmupSteps
   Number of MC steps discarded as a warmup or equilibration period of
   the random walk. If this is too small, it will bias the optimization
   procedure.

blocks
   Number of average energy values written to output files. Should be
   greater than 200 for meaningful statistical analysis of output data
   (e.g., via ``qmca``).

subSteps
   Number of MC steps in between energy evaluations. Each energy
   evaluation is expensive, so taking a few steps to decorrelate between
   measurements can be more efficient. Will be less efficient with many
   substeps.

nonlocalpp,useBuffer
   If ``nonlocalpp="no,"`` then the nonlocal part of the pseudopotential
   is not included when computing the cost function. If
   ``useBuffer="yes,"`` then temporary data is stored to speed up
   nonlocal pseudopotential evaluation at the expense of memory
   consumption.

loop max
   Number of times to repeat the optimization. Using the resulting
   wavefunction from the previous optimization in the next one improves
   the results. Typical choices range between 8 and 16.

The cost function defines the quantity to be minimized during
optimization. The three components of the cost function, energy,
unreweighted variance, and reweighted variance should sum to one.
Dedicating 100% of the cost function to unreweighted variance is often a
good choice. Another common choice is to try 90/10 or 80/20 mixtures of
reweighted variance and energy. Using 100% energy minimization is
desirable for reducing DMC pseudopotential localization errors, but the
optimization process is less stable and should be attempted only after
first performing several cycles of, for example, variance minimization
(the entire ``loop`` section can be duplicated with a different cost
function each time).

Replace ``MAX``, ``EVCOST``, ``UVCOST``, ``RVCOST``, ``TS``, and ``SAMPLES`` in the ``loop`` with appropriate starting values in the ``O.q0.opt.in.xml`` input file.  Perform the optimization run by typing

::

  mpirun -np 4 qmcpack O.q0.opt.in.xml >&O.q0.opt.out&

.. jobrun_vesta qmcpack O.q0.opt.in.xml

The run should take only a few minutes for reasonable values of loop ``max`` and ``samples``.

The log file output will appear in ``O.q0.opt.out``.  The beginning of each linear optimization will be marked with text similar to

::

  =========================================================
    Start QMCFixedSampleLinearOptimize
    File Root O.q0.opt.s011 append = no
  =========================================================

At the end of each optimization section the change in cost function, new values for the Jastrow parameters, and elapsed wall clock time are reported:

::

  OldCost: 7.0598901869e-01 NewCost: 7.0592576381e-01 Delta Cost:-6.3254886314e-05
  ...
   <optVariables href="O.q0.opt.s011.opt.xml">
  uu_0 6.9392504232e-01 1 1  ON 0
  uu_1 4.9690781460e-01 1 1  ON 1
  uu_2 4.0934542375e-01 1 1  ON 2
  uu_3 3.7875640157e-01 1 1  ON 3
  uu_4 3.7308380014e-01 1 1  ON 4
  uu_5 3.5419786809e-01 1 1  ON 5
  uu_6 4.3139019377e-01 1 1  ON 6
  uu_7 1.9344371667e-01 1 1  ON 7
  ud_0 3.9219009713e-01 1 1  ON 8
  ud_1 1.2352664647e-01 1 1  ON 9
  ud_2 4.4048945133e-02 1 1  ON 10
  ud_3 2.1415676741e-02 1 1  ON 11
  ud_4 1.5201803731e-02 1 1  ON 12
  ud_5 2.3708169445e-02 1 1  ON 13
  ud_6 3.4279064930e-02 1 1  ON 14
  ud_7 4.3334583596e-02 1 1  ON 15
  eO_0 -7.8490123937e-01 1 1  ON 16
  eO_1 -6.6726618338e-01 1 1  ON 17
  eO_2 -4.8753453838e-01 1 1  ON 18
  eO_3 -3.0913993774e-01 1 1  ON 19
  eO_4 -1.7901872177e-01 1 1  ON 20
  eO_5 -8.6199000697e-02 1 1  ON 21
  eO_6 -4.0601160841e-02 1 1  ON 22
  eO_7 -4.1358075061e-03 1 1  ON 23
   </optVariables>
  ...
   QMC Execution time = 2.8218972974e+01 secs

The cost function should decrease during each linear optimization (``Delta cost < 0``).  Try "``grep OldCost *opt.out.``"  You should see something like this:

::

  OldCost: 1.2655186572e+00 NewCost: 7.2443875597e-01 Delta Cost:-5.4107990118e-01
  OldCost: 7.2229830632e-01 NewCost: 6.9833678217e-01 Delta Cost:-2.3961524143e-02
  OldCost: 8.0649629434e-01 NewCost: 8.0551871147e-01 Delta Cost:-9.7758287036e-04
  OldCost: 6.6821241388e-01 NewCost: 6.6797703487e-01 Delta Cost:-2.3537901148e-04
  OldCost: 7.0106275099e-01 NewCost: 7.0078055426e-01 Delta Cost:-2.8219672877e-04
  OldCost: 6.9538522411e-01 NewCost: 6.9419186712e-01 Delta Cost:-1.1933569922e-03
  OldCost: 6.7709626744e-01 NewCost: 6.7501251165e-01 Delta Cost:-2.0837557922e-03
  OldCost: 6.6659923822e-01 NewCost: 6.6651737755e-01 Delta Cost:-8.1860671682e-05
  OldCost: 7.7828995609e-01 NewCost: 7.7735482525e-01 Delta Cost:-9.3513083900e-04
  OldCost: 7.2717974404e-01 NewCost: 7.2715201115e-01 Delta Cost:-2.7732880747e-05
  OldCost: 6.9400639873e-01 NewCost: 6.9257183689e-01 Delta Cost:-1.4345618444e-03
  OldCost: 7.0598901869e-01 NewCost: 7.0592576381e-01 Delta Cost:-6.3254886314e-05

Blocked averages of energy data, including the kinetic energy and components of the potential energy, are written to ``scalar.dat`` files.  The first is named "``O.q0.opt.s000.scalar.dat``," with a series number of zero (``s000``).  In the end there will be ``MAX`` of them, one for each series.

When the job has finished, use the ``qmca`` tool to assess the effectiveness of the optimization process.  To look at just the total energy and the variance, type "``qmca -q ev O.q0.opt*scalar*``."  This will print the energy, variance, and the variance/energy ratio in Hartree units:

::

                    LocalEnergy               Variance           ratio
  O.q0.opt  series 0  -15.739585 +/- 0.007656   0.887412 +/- 0.010728   0.0564
  O.q0.opt  series 1  -15.848347 +/- 0.004089   0.318490 +/- 0.006404   0.0201
  O.q0.opt  series 2  -15.867494 +/- 0.004831   0.292309 +/- 0.007786   0.0184
  O.q0.opt  series 3  -15.871508 +/- 0.003025   0.275364 +/- 0.006045   0.0173
  O.q0.opt  series 4  -15.865512 +/- 0.002997   0.278056 +/- 0.006523   0.0175
  O.q0.opt  series 5  -15.864967 +/- 0.002733   0.278065 +/- 0.004413   0.0175
  O.q0.opt  series 6  -15.869644 +/- 0.002949   0.273497 +/- 0.006141   0.0172
  O.q0.opt  series 7  -15.868397 +/- 0.003838   0.285451 +/- 0.007570   0.0180
  ...

Plots of the data can also be obtained with the “``-p``” option
(“``qmca -p -q ev O.q0.opt*scalar*``”).

Identify which optimization series is the "best" according to your cost function.  It is likely that multiple series are similar in quality.  Note the ``opt.xml`` file corresponding to this series.  This file contains the final value of the optimized Jastrow parameters to be used in the DMC calculations of the next section of the lab.

**Questions and Exercises**

#. What is the acceptance ratio of your optimization runs? (use
   “texttqmca -q ar O.q0.opt*scalar\*”) Do you expect the MC sampling to
   be efficient?

#. How do you know when the optimization process has converged?

#. (optional) Optimization is sometimes sensitive to initial guesses of
   the parameters. If you have time, try varying the initial parameters,
   including the cutoff radius (``rcut``) of the Jastrow factors
   (remember to change ``id`` in the ``<project/>`` element). Do you
   arrive at a similar set of final Jastrow parameters? What is the
   lowest variance you are able to achieve?

DMC timestep extrapolation I: neutral oxygen atom
-------------------------------------------------

The DMC algorithm contains two biases in addition to the fixed node and pseudopotential approximations that are important to control: time step and population control bias.  In this section we focus on estimating and removing time step bias from DMC calculations.  The essential fact to remember is that the bias vanishes as the time step goes to zero, while the needed computer time increases inversely with the time step.

In the same directory you used to perform wavefunction optimization (``oxygen_atom``) you will find a sample DMC input file for the neutral oxygen atom named ``O.q0.dmc.in.xml``.  Open this file in a text editor and note the differences from the optimization case.  Wavefunction information is no longer included from ``pw2qmcpack`` but instead should come from the optimization run:

::

  <!-- OPT_XML is from optimization, e.g. O.q0.opt.s008.opt.xml -->
  <include href="OPT_XML"/>

Replace "``OPT_XML``" with the ``opt.xml`` file corresponding to the best Jastrow parameters you found in the last section (this is a file name similar to ``O.q0.opt.s008.opt.xml``).

The QMC calculation section at the bottom is also different.  The linear optimization blocks have been replaced with XML describing a VMC run followed by DMC.  Descriptions of the input keywords follow.

timestep
   Time step of the VMC/DMC random walk. In VMC choose a time step
   corresponding to an acceptance ratio of about 50%. In DMC the
   acceptance ratio is often above 99%.

warmupSteps
   Number of MC steps discarded as a warmup or equilibration period of
   the random walk.

steps
   Number of MC steps per block. Physical quantities, such as the total
   energy, are averaged over walkers and steps.

blocks
   Number of blocks. This is also the number of average energy values
   written to output files. The number should be greater than 200 for
   meaningful statistical analysis of output data (e.g., via ``qmca``).
   The total number of MC steps each walker takes is
   ``blocks``\ :math:`\times`\ ``steps``.

samples
   VMC only. This is the number of walkers used in subsequent DMC runs.
   Each DMC walker is initialized with electron positions sampled from
   the VMC random walk.

nonlocalmoves
   DMC only. If yes/no, use the locality approximation/T-moves for
   nonlocal pseudopotentials. T-moves generally improve the stability of
   the algorithm and restore the variational principle for small systems
   (T-moves version 1).

The purpose of the VMC run is to provide initial electron positions for each DMC walker.  Setting :math:`\texttt{walkers}=1` in the VMC block ensures there will be only one VMC walker per execution thread.  There will be a total of 4 VMC walkers in this case (see ``O.q0.dmc.qsub.in``).  We want the electron positions used to initialize the DMC walkers to be decorrelated from one another.  A VMC walker will often decorrelate from its current position after propagating for a few Ha :math:`^{-1}` in imaginary time (in general, this is system dependent).  This leads to a rough rule of thumb for choosing ``blocks`` and ``steps`` for the VMC run (:math:`\texttt{vwalkers}=4` here):

.. math::
  :label: eq68

  \begin{aligned}
     \texttt{VBLOCKS}\times\texttt{VSTEPS} \ge \frac{\texttt{DWALKERS}}{\texttt{VWALKERS}} \frac{5~\textrm{Ha}^{-1}}{\texttt{VTIMESTEP}}\end{aligned}

Fill in the VMC XML block with appropriate values for these parameters.
There should be more than one DMC walker per thread and enough walkers
in total to avoid population control bias. The general rule of thumb is
to have more than :math:`\sim 2,000` walkers, although the dependence of
the total energy on population size should be explicitly checked from
time to time.

To study time step bias, we will perform a sequence of DMC runs over a
range of time steps (:math:`0.1` Ha\ :math:`^{-1}` is too large, and
time steps below :math:`0.002` Ha\ :math:`^{-1}` are probably too
small). A common approach is to select a fairly large time step to begin
with and then decrease the time step by a factor of two in each
subsequent DMC run. The total amount of imaginary time the walker
population propagates should be the same for each run. A simple way to
accomplish this is to choose input parameters in the following way

.. math::
  :label: eq69


   \begin{aligned}
     \texttt{timestep}_{n}    &= \texttt{timestep}_{n-1}/2\nonumber\\
     \texttt{warmupSteps}_{n} &= \texttt{warmupSteps}_{n-1}\times 2\nonumber\\
     \texttt{blocks}_{n}      &= \texttt{blocks}_{n-1}\nonumber\\
     \texttt{steps}_{n}       &= \texttt{steps}_{n-1}\times 2\end{aligned}
