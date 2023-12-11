.. _lab-condensed-matter:

Lab 4: Condensed Matter Calculations
====================================



.. **Lab author: Luke Shulenburger**
   Sandia National Laboratories is a multiprogram laboratory managed and
   operated by Sandia Corporation, a wholly owned subsidiary of Lockheed
   Martin Corporation, for the U.S. Department of Energy’s National
   Nuclear Security Administration under Contract No. DE-AC04-94AL85000.

Topics covered in this lab
--------------------------

-  Tiling DFT primitive cells into QMC supercells

-  Reducing finite-size errors via extrapolation

-  Reducing finite-size errors via averaging over twisted boundary
   conditions

-  Using the B-spline mesh factor to reduce memory requirements

-  Using a coarsely resolved vacuum buffer region to reduce memory
   requirements

-  Calculating the DMC total energies of representative 2D and 3D
   extended systems

Lab directories and files
-------------------------

::

  labs/lab4_condensed_matter/
  ├── Be-2at-setup.py           - DFT only for prim to conv cell
  ├── Be-2at-qmc.py             - QMC only for prim to conv cell
  ├── Be-16at-qmc.py            - DFT and QMC for prim to 16 atom cell
  ├── graphene-setup.py         - DFT and OPT for graphene
  ├── graphene-loop-mesh.py     - VMC scan over orbital bspline mesh factors
  ├── graphene-final.py         - DMC for final meshfactor
  └── pseudopotentials          - pseudopotential directory
      ├── Be.ncpp                 - Be PP for Quantum ESPRESSO
      ├── Be.xml                  - Be PP for QMCPACK
      ├── C.BFD.upf               - C  PP for Quantum ESPRESSO
      └── C.BFD.xml               - C  PP for QMCPACK

The goal of this lab is to introduce you to the somewhat specialized problems involved in performing DMC calculations on condensed matter as opposed to the atoms and molecules that were the focus of the preceding labs.   Calculations will be performed on two different systems.  Firstly, we will perform a series of calculations on BCC beryllium, focusing on the necessary methodology to limit finite-size effects.  Secondly, we will perform calculations on graphene as an example of a system where QMCPACK’s capability to handle cases with mixed periodic and open boundary conditions is useful.  This example will also focus on strategies to limit memory usage for such systems.
All of the calculations performed in this lab will use the Nexus workflow management system, which vastly simplifies the process by automating the steps of generating trial wavefunctions and performing DMC calculations.

Preliminaries
-------------

For any DMC calculation, we must start with a trial wavefunction. As is typical for our calculations of condensed matter, we will produce this wavefunction using DFT.  Specifically, we will use QE to generate a Slater determinant of SPOs.  This is done as a three-step process.  First, we calculate the converged charge density by performing a DFT calculation with a fine grid of k-points to fully sample the Brillouin zone.  Next, a non-self- consistent calculation is performed at the specific k-points needed for the supercell and twists needed in the DMC calculation (more on this later).  Finally, a wavefunction is converted from the binary representation used by QE to the portable hdf5 representation used by QMCPACK.

The choice of k-points necessary to generate the wavefunctions depends on both the supercell chosen for the DMC calculation and by the supercell twist vectors needed.  Recall that the wavefunction in a plane-wave DFT calculation is written using Bloch's theorem as:

.. math::
  :label: eq70

  \Psi(\vec{r}) = e^{i\vec{k}\cdot\vec{r}}u(\vec{r})\:,

where :math:`\vec{k}` is confined to the first Brillouin zone of the
cell chosen and :math:`u(\vec{r})` is periodic in this simulation cell.
A plane-wave DFT calculation stores the periodic part of the
wavefunction as a linear combination of plane waves for each SPO at all
k-points selected. The symmetry of the system allows us to generate an
arbitrary supercell of the primitive cell as follows: Consider the set
of primitive lattice vectors, :math:`\{ \mathbf{a}^p_1, \mathbf{a}^p_2,
\mathbf{a}^p_3\}`. We may write these vectors in a matrix,
:math:`\mathbf{L}_p`, the rows of which are the primitive lattice
vectors. Consider a nonsingular matrix of integers, :math:`\mathbf{S}`.
A corresponding set of supercell lattice vectors,
:math:`\{\mathbf{a}^s_1, \mathbf{a}^s_2, \mathbf{a}^s_3\}`, can be
constructed by the matrix product

.. math::
  :label: eq71

  \mathbf{a}^s_i = S_{ij} \mathbf{a}^p_j]\:.

If the primitive cell contains :math:`N_p` atoms, the supercell will
then contain :math:`N_s = |\det(\mathbf{S})| N_p` atoms.

Now, the wavefunciton at any point in this new supercell can be related
to the wavefunction in the primitive cell by finding the linear
combination of primitive lattice vectors that maps this point back to
the primitive cell:

.. math::
  :label: eq72

  \vec{r}' = \vec{r} + x \mathbf{a}^p_1 + y \mathbf{a}^p_2 + z\mathbf{a}^p_3 = \vec{r} + \vec{T}\:,

where :math:`x, y, z` are integers. Now the wavefunction in the
supercell at point :math:`\vec{r}'` can be written in terms of the
wavefunction in the primitive cell at :math:`\vec{r}'` as:

.. math::
  :label: eq73

  \Psi(\vec{r}) = \Psi(\vec{r}') e^{i \vec{T} \cdot \vec{k}}\:,


where :math:`\vec{k}` is confined to the first Brillouin zone of the
primitive cell. We have also chosen the supercell twist vector, which
places a constraint on the form of the wavefunction in the supercell.
The combination of these two constraints allows us to identify family of
N k-points in the primitive cell that satisfy the constraints. Thus, for
a given supercell tiling matrix and twist angle, we can write the
wavefunction everywhere in the supercell by knowing the wavefunction a N
k-points in the primitive cell. This means that the memory necessary to
store the wavefunction in a supercell is only linear in the size of the
supercell rather than the quadratic cost if symmetry were neglected.

Total energy of BCC beryllium
-----------------------------

When performing calculations of periodic solids with QMC, it is essential to work with a reasonable size supercell rather than the primitive cells that are common in mean field calculations.  Specifically, all of the finite-size correction schemes discussed in the morning require that the exchange-correlation hole be considerably smaller than the periodic simulation cell.  Additionally, finite-size effects are lessened as the distance between the electrons in the cell and their periodic images increases, so it is advantageous to generate supercells that are as spherical as possible to maximize this distance.  However, a competing consideration is that when calculating total energies we often want to extrapolate the energy per particle to the thermodynamic limit by means of the following formula in three dimensions:

.. math::
  :label: eq74

   E_{\inf} = C + E_{N}/N\:.

This formula derived assuming the shape of the supercells is consistent
(more specifically that the periodic distances scale uniformly with
system size), meaning we will need to do a uniform tiling, that is,
:math:`2\times2\times2`, :math:`3\times3\times3`, etc. As a
:math:`3\times3\times3` tiling is 27 times larger than the supercell and
the practical limit of DMC is on the order of 200 atoms (depending on
Z), sometimes it is advantageous to choose a less spherical supercell
with fewer atoms rather than a more spherical one that is too expensive
to tile.

In the case of a BCC crystal, it is possible to tile the one atom
primitive cell to a cubic supercell only by doubling the number of
electrons. This is the best possible combination of a small number of
atoms that can be tiled and a regular box that maximizes the distance
between periodic images. We will need to determine the tiling matrix S
that generates this cubic supercell by solving the following equation
for the coefficients of the S matrix:

.. math::
  :label: eq75

  \left[\begin{array}{rrr}
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 1
    \end{array}\right] =  \left[\begin{array}{rrr}
    s_{11} & s_{12} & s_{13} \\
    s_{21} & s_{22} & s_{23} \\
    s_{31} & s_{32} & s_{33}
    \end{array}\right] \cdot
  \left[\begin{array}{rrr}
    0.5 &  0.5 & -0.5 \\
   -0.5 &  0.5 &  0.5 \\
    0.5 & -0.5 &  0.5
  \end{array}\right]\:.

We will now use Nexus to generate the trial wavefunction for this BCC beryllium.

Fortunately, the Nexus will handle determination of the proper k-vectors given the tiling matrix.  All that is needed is to place the tiling matrix in the ``Be-2at-setup.py`` file.   Now the definition of the physical system is

::

  bcc_Be = generate_physical_system(
      lattice    = 'cubic',
      cell       = 'primitive',
      centering  = 'I',
      atoms      = 'Be',
      constants  = 3.490,
      units      = 'A',
      net_charge = 0,
      net_spin   = 0,
      Be         = 2,
      tiling     = [[a,b,c],[d,e,f],[g,h,i]],
      kgrid      = kgrid,
      kshift     = (.5,.5,.5)
      )

where the tiling line should be replaced with the preceding row major
tiling matrix. This script file will now perform a converged DFT
calculation to generate the charge density in a directory called
``bcc-beryllium/scf`` and perform a non-self-consistend DFT calculation
to generate SPOs in the directory ``bcc-beryllium/nscf``. Fortunately,
Nexus will calculate the required k-points needed to tile the
wavefunction to the supercell, so all that is necessary is the
granularity of the supercell twists and whether this grid is shifted
from the origin. Once this is finished, it performs the conversion from
pwscf’s binary format to the hdf5 format used by QMCPACK. Finally, it
will optimize the coefficients of 1-body and 2-body Jastrow factors in
the supercell defined by the tiling matrix.

Run these calculations by executing the script ``Be-2at-setup.py``. You
will notice the small calculations required to generate the wavefunction
of beryllium in a one-atom cell are rather inefficient to run on a
high-performance computer such as vesta in terms of the time spent doing
calculations versus time waiting on the scheduler and booting compute
nodes. One of the benefits of the portable HDF format that is used by
QMCPACK is that you can generate data like wavefunctions on a local
workstation or other convenient resource and use high-performance
clusters for the more expensive QMC calculations.

In this case, the wavefunction is generated in the directory
``bcc-beryllium/nscf-2at_222/pwscf_ output`` in a file called
``pwscf.pwscf.h5``. For debugging purposes, it can be useful to verify
that the contents of this file are what you expect. For instance, you
can use the tool ``h5ls`` to check the geometry of the cell where the
DFT calculations were performed or the number of k-points or electrons
in the calculation. This is done with the command h5ls -d
pwscf.pwscf.h5/supercell or h5ls -d pwscf.pwscf.h5/electrons.

In the course of running ``Be-2at-setup.py``, you will get an error when
attempting to perform the VMC and wavefunction optimization
calculations. This is because the wavefunction has generated supercell
twists of the form (+/- 1/4, +/- 1/4, +/- 1/4). In the case that the
supercell twist contains only 0 or 1/2, it is possible to operate
entirely with real arithmetic. The executable that has been indicated in
``Be-2at-setup.py`` was compiled for this case. Note that where
possible, the memory use is a factor of two less than the general case
and the calculations are somewhat faster. However, it is often necessary
to perform calculations away from these special twist angles to reduce
finite-size effects. To fix this, delete the directory
``bcc-beryllium/opt-2at``, change the line near the top of
``Be-2at-setup.py`` from

::

  qmcpack    = '/soft/applications/qmcpack/Binaries/qmcpack'

to

::

  qmcpack    = '/soft/applications/qmcpack/Binaries/qmcpack_comp'

and rerun the script.

When the optimization calculation has finished, check that everything has proceeded correctly by looking at the output in the ``opt-2at`` directory.  Firstly, you can grep the output file for Delta to see if the cost function has indeed been decreasing during the optimization.  You should find something like this:

::

  OldCost: 4.8789147e-02 NewCost: 4.0695360e-02 Delta Cost:-8.0937871e-03
  OldCost: 3.8507795e-02 NewCost: 3.8338486e-02 Delta Cost:-1.6930674e-04
  OldCost: 4.1079105e-02 NewCost: 4.0898345e-02 Delta Cost:-1.8076319e-04
  OldCost: 4.2681333e-02 NewCost: 4.2356598e-02 Delta Cost:-3.2473514e-04
  OldCost: 3.9168577e-02 NewCost: 3.8552883e-02 Delta Cost:-6.1569350e-04
  OldCost: 4.2176276e-02 NewCost: 4.2083371e-02 Delta Cost:-9.2903058e-05
  OldCost: 4.3977361e-02 NewCost: 4.2865751e-02 Delta Cost:-1.11161830-03
  OldCost: 4.1420944e-02 NewCost: 4.0779569e-02 Delta Cost:-6.4137501e-04

which shows that the starting wavefunction was fairly good and that most
of the optimization occurred in the first step. Confirm this by using
``qmca`` to look at how the energy and variance changed over the course
of the calculation with the command: ``qmca -q ev -e 10 *.scalar.dat``
executed in the ``opt-2at directory``. You should get output like the
following:

::

                    LocalEnergy               Variance             ratio
  opt  series 0  -2.159139 +/- 0.001897   0.047343 +/- 0.000758   0.0219
  opt  series 1  -2.163752 +/- 0.001305   0.039389 +/- 0.000666   0.0182
  opt  series 2  -2.160913 +/- 0.001347   0.040879 +/- 0.000682   0.0189
  opt  series 3  -2.162043 +/- 0.001223   0.041183 +/- 0.001250   0.0190
  opt  series 4  -2.162441 +/- 0.000865   0.039597 +/- 0.000342   0.0183
  opt  series 5  -2.161287 +/- 0.000732   0.039954 +/- 0.000498   0.0185
  opt  series 6  -2.163458 +/- 0.000973   0.044431 +/- 0.003583   0.0205
  opt  series 7  -2.163495 +/- 0.001027   0.040783 +/- 0.000413   0.0189

Now that the optimization has completed successfully, we can perform DMC
calculations. The first goal of the calculations will be to try to
eliminate the 1-body finite-size effects by twist averaging. The script
``Be-2at-qmc.py`` has the necessary input. Note that on line 42 two
twist grids are specified, (2,2,2) and (3,3,3). Change the tiling matrix
in this input file as in ``Be-2at-qmc.py`` and start the calculations.
Note that this workflow takes advantage of QMCPACK’s capability to group
jobs. If you look in the directory ``dmc-2at_222`` at the job submission
script (``dmc.qsub.in``), you will note that rather than operating on an
XML input file, ``qmcapp`` is targeting a text file called ``dmc.in``.
This file is a simple text file that contains the names of the eight XML
input files needed for this job, one for each twist. When operated in
this mode, QMCPACK will use MPI groups to run multiple copies of itself
within the same MPI context. This is often useful both in terms of
organizing calculations and for taking advantage of the large job sizes
that computer centers often encourage.

The DMC calculations in this case are designed to complete in a few
minutes. When they have finished running, first look at the
``scalar.dat`` files corresponding to the DMC calculations at the
various twists in ``dmc-2at_222``. Using a command such as
``qmca -q ev -e 32 *.s001.scalar.dat`` (with a suitably chosen number of
blocks for the equilibration), you will see that the DMC energy in each
calculation is nearly identical within the statistical uncertainty of
the calculations. In the case of a large supercell, this is often
indicative of a situation where the Brillouin zone is so small that the
1-body finite-size effects are nearly converged without any twist
averaging. In this case, however, this is because of the symmetry of the
system. For this cubic supercell, all of the twist angles chosen in this
shifted :math:`2\times2\times2` grid are equivalent by symmetry. In the
case where substantial resources are required to equilibrate the DMC
calculations, it can be beneficial to avoid repeating such twists and
instead simply weight them properly. In this case, however, where the
equilibration is inexpensive, there is no benefit to adding such
complexity as the calculations can simply be averaged together and the
result is equivalent to performing a single longer calculation.

Using the command ``qmc -a -q ev -e 16 *.s001.scalar.dat``, average the
DMC energies in ``dmc-2at_222 and dmc-2at_333`` to see whether the
1-body finite-size effects are converged with a :math:`3\times3\times3`
grid of twists. When using beryllium as a metal, the convergence is
quite poor (0.025 Ha/Be or 0.7 eV/Be). If this were a production
calculation it would be necessary to perform calculations on much larger
grids of supercell twists to eliminate the 1-body finite-size effects.

In this case there are several other calculations that would warrant a
high priority. Script ``Be-16at-qmc.py`` has been provided in which you
can input the appropriate tiling matrix for a 16-atom cell and perform
calculations to estimate the 2-body finite-size effects, which will also
be quite large in the 2-atom calculations. This script will take
approximately 30 minutes to run to completion, so depending on your
interest, you can either run it or work to modify the scripts to address
the other technical issues that would be necessary for a production
calculation such as calculating the population bias or the time step
error in the DMC calculations.

Another useful exercise would be to attempt to validate this PP by
calculating the ionization potential and electron affinity of the
isolated atom and compare it with the experimental values: IP = 9.3227
eV , EA = 2.4 eV.

Handling a 2D system: graphene
------------------------------

In this section we examine a calculation of an isolated sheet of
graphene. Because graphene is a 2D system, we will take advantage of
QMCPACK’s capability to mix periodic and open boundary conditions to
eliminate and spurious interaction of the sheet with its images in the z
direction. Run the script ``graphene-setup.py``, which will generate the
wavefunction and optimize one and two body jastrow factors. In the
script; notice line 160: bconds = ’ppn’ in the generate_qmcpack
function, which specifies this mix of open and periodic boundary
conditions. Consequently, the atoms will need to be kept away from this
open boundary in the z direction as the electronic wavefunction will not
be defined outside of the simulation box in this direction. For this
reason, all of the atom positions at the beginning of the file have z
coordinates 7.5. At this point, run the script ``graphene-setup.py``.

Aside from the change in boundary conditions, the main thing that
distinguishes this kind of calculation from the previous beryllium
example is the large amount of vacuum in the cell. Although this is a
very small calculation designed to run quickly in the tutorial, in
general a more converged calculation would quickly become memory limited
on an architecture like BG/Q. When the initial wavefunction optimization
has completed to your satisfaction, run the script
``graphene-loop-mesh.py``. This examines within VMC an approach to
reducing the memory required to store the wavefunction. In
``graphene-loop-mesh.py``, the spacing between the B-spline points is
varied uniformly. The mesh spacing is a prefactor to the linear spacing
between the spline points, so the memory use goes as the cube of the
meshfactor. When you run the calculations, examine the
``.s000.scalar.dat`` files with ``qmca`` to determine the lowest
possible mesh spacing that preserves both the VMC energy and the
variance.

Finally, edit the file ``graphene-final.py``, which will perform two DMC
calculations. In the first, (qmc1) replace the following lines:

::

  meshfactor   = xxx,
  precision    = '---',

with the values you have determined will perform the calculation with as small as possible wavefunction.  Note that we can also use single precision arithmetic to store the wavefunction by specifying precision=`single.'  When you run the script, compare the output of the two DMC calculations in terms of energy and variance.  Also, see if you can calculate the fraction of memory that you were able to save by using a meshfactor other than 1 and single precision arithmetic.

Conclusion
----------

Upon completion of this lab, you should be able to use Nexus to perform DMC calculations on periodic solids when provided with a PP.  You should also be able to reduce the size of the wavefunction in a solid-state calculation in cases where memory is a limiting factor.
