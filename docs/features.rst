.. _chap:features:

Features of QMCPACK
===================

Production-level features
-------------------------

The following list contains the main production-level features of
QMCPACK. If you do not see a specific feature that you are interested
in, see the remainder of this manual and ask whether that specific
feature is available or can be quickly brought to the full production
level.

-  Variational Monte Carlo (VMC)

-  Diffusion Monte Carlo (DMC)

-  Reptation Monte Carlo

-  Single and multideterminant Slater Jastrow wavefunctions

-  Wavefunction updates using optimized multideterminant algorithm of
   Clark et al.

-  Backflow wavefunctions

-  One, two, and three-body Jastrow factors

-  Excited state calculations via flexible occupancy assignment of
   Slater determinants

-  All electron and nonlocal pseudopotential calculations

-  Casula T-moves for variational evaluation of nonlocal
   pseudopotentials (non-size-consistent and size-consistent variants)

-  Wavefunction optimization using the “linear method” of Umrigar and
   coworkers, with arbitrary mix of variance and energy in the objective
   function

-  Blocked, low memory adaptive shift optimizer of Zhao and Neuscamman

-  Gaussian, Slater, plane-wave, and real-space spline basis sets for
   orbitals

-  Interface and conversion utilities for plane-wave wavefunctions from
   Quantum Espresso (Plane-Wave Self-Consistent Field package [PWSCF])

-  Interface and conversion utilities for Gaussian-basis wavefunctions
   from GAMESS

-  Easy extension and interfacing to other electronic structure codes
   via standardized XML and HDF5 inputs

-  MPI parallelism

-  Fully threaded using OpenMP

-  GPU (NVIDIA CUDA) implementation (limited functionality)

-  HDF5 input/output for large data

-  Nexus: advanced workflow tool to automate all aspects of QMC
   calculation from initial DFT calculations through to final analysis

-  Analysis tools for minimal environments (Perl only) through to
   Python-based environments with graphs produced via matplotlib
   (included with Nexus)

SoA optimizations and improved algorithms
-----------------------------------------

The Structure-of-Arrays (SoA) optimizations
:cite:`IPCC_SC17` are a set of improved data layouts
facilitating vectorization on modern CPUs with wide SIMD units. **For
many calculations and architectures, the SoA implementation more than
doubles the speed of the code.** This so-called SoA implementation
replaces the older, less efficient Array-of-Structures (AoS) code and
can be enabled or disabled at compile time. The memory footprint is also
reduced in the SoA implementation by better algorithms, enabling more
systems to be run.

The SoA build was made the default for QMCPACK v3.7.0. As described in :ref:`cmakeoptions`, the SoA
implementation can be disabled by configuring with ``-DENABLE_SOA=0``.

The SoA code path currently does *not* support:

-  Backflow wavefunctions

-  Many observables

The code should abort with a message referring to AoS vs SoA features if
any unsupported feature is invoked. In this case the AoS build should be
used by configuring with ``-DENABLE_SOA=0``. In addition, please inform the developers via
GitHub or Google Groups so that porting these features can be
prioritized.

Core features are heavily tested in both SoA and AoS versions. If using
untested and noncore features in the SoA code, please compare the AoS
and SoA results carefully.

Supported GPU features
----------------------

The GPU implementation supports multiple GPUs per node, with MPI tasks
assigned in a round-robin order to the GPUs. Currently, for large runs,
1 MPI task should be used per GPU per node. For smaller calculations,
use of multiple MPI tasks per GPU might yield improved performance.
Supported GPU features:

-  VMC, wavefunction optimization, DMC.

-  Periodic and open boundary conditions. Mixed boundary conditions are
   not yet supported.

-  Wavefunctions:

   #. Single Slater determinants with 3D B-spline orbitals.
      Twist-averaged boundary conditions and complex wavefunctions are
      fully supported. Gaussian type orbitals are not yet supported.

   #. Hybrid mixed basis representation in which orbitals are
      represented as 1D splines times spherical harmonics in spherical
      regions (muffin tins) around atoms and 3D B-splines in the
      interstitial region.

   #. One-body and two-body Jastrows represented as 1D B-splines.
      Three-body Jastrow functions are not yet supported.

-  Semilocal (nonlocal and local) pseudopotentials, Coulomb interaction
   (electron-electron, electron-ion), and model periodic Coulomb (MPC)
   interaction.

Beta test features
------------------

This section describes developmental features in QMCPACK that might be
ready for production but that require additional testing, features, or
documentation to be ready for general use. We describe them here because
they offer significant benefits and are well tested in specific cases.

Auxiliary-Field Quantum Monte Carlo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The orbital-space Auxiliary-Field Quantum Monte Carlo (AFMQC) method is
now available in QMCPACK. The main input for the code is the matrix
elements of the Hamiltonian in a given single particle basis set, which
must be produced from mean-field calculations such as Hartree-Fock or
density functional theory. The code and many features are in
development. Check the latest version of QMCPACK for an up-to-date
description of available features. A partial list of the current
capabilities of the code follows. For a detailed description of the
available features, see  :ref:`afqmc`.

-  Phaseless AFQMC algorithm of Zhang et al. (S. Zhang and H. Krakauer.
   2003. “Quantum Monte Carlo Method using Phase-Free Random Walks with
   Slater Determinants." *PRL* 90: 136401).

-  “Hybrid" and “local energy" propagation schemes.

-  Hamiltonian matrix elements from (1) Molpro’s FCIDUMP format (which
   can be produced by Molpro, PySCF, and VASP) and (2) internal HDF5
   format produced by PySCF (see AFQMC section below).

-  AFQMC calculations with RHF (closed-shell doubly occupied), ROHF
   (open-shell doubly occupied), and UHF (spin polarized broken
   symmetry) symmetry.

-  Single and multideterminant trial wavefunctions. Multideterminant
   expansions with either orthogonal or nonorthogonal determinants.

-  Fast update scheme for orthogonal multideterminant expansions.

-  Distributed propagation algorithms for large systems. Enables
   calculations where data structures do not fit on a single node.

-  Complex implementation for PBC calculations with complex integrals.

-  Sparse representation of large matrices for reduced memory usage.

-  Mixed and back-propagated estimators.

-  Specialized implementation for solids with k-point symmetry (e.g.
   primitive unit cells with kpoints).

-  Efficient GPU implementation (currently limited to solids with
   k-point symmetry).

Sharing of spline data across multiple GPUs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sharing of GPU spline data enables distribution of the data across
multiple GPUs on a given computational node. For example, on a
two-GPU-per-node system, each GPU would have half of the orbitals. This
allows use of larger overall spline tables than would fit in the memory
of individual GPUs and potentially up to the total GPU memory on a node.
To obtain high performance, large electron counts or a high-performing
CPU-GPU interconnect is required.

To use this feature, the following needs to be done:

-  The CUDA Multi-Process Service (MPS) needs to be used (e.g., on OLCF
   Summit/SummitDev use “-alloc_flags gpumps" for bsub). If MPI is not
   detected, sharing will be disabled.

-  CUDA_VISIBLE_DEVICES needs to be properly set to control each rank’s
   visible CUDA devices (e.g., on OLCF Summit/SummitDev create a
   resource set containing all GPUs with the respective number of ranks
   with “jsrun –task-per-rs Ngpus -g Ngpus").

-  In the determinant set definition of the <wavefunction> section, the
   “gpusharing" parameter needs to be set (i.e., <determinantset
   gpusharing=“yes">). See
   :ref:`spo-spline`.

.. bibliography:: bibliography.bib
