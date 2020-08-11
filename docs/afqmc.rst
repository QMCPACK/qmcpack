.. _afqmc:

Auxiliary-Field Quantum Monte Carlo
===================================

The AFQMC method is an orbital-space formulation of the imaginary-time propagation algorithm. We refer the reader to one of the review articles on the method :cite:`AFQMC_review,PhysRevLett.90.136401,PhysRevE.70.056702` for a detailed description of the algorithm. It uses the Hubbard-Stratonovich transformation to express the imaginary-time propagator, which is inherently a 2-body operator, as an integral over 1-body propagators, which can be efficiently applied to an arbitrary Slater determinant. This transformation allows us to represent the interacting many-body system as an average over a noninteracting system (e.g., Slater determinants) in a time-dependent fluctuating external field (the Auxiliary fields). The walkers in this case represent nonorthogonal Slater determinants, whose time average represents the desired quantum state. QMCPACK currently implements the phaseless AFQMC algorithm of Zhang and Krakauer :cite:`PhysRevLett.90.136401`, where a trial wavefunction is used to project the simulation to the real axis, controlling the fermionic sign problem at the expense of a bias. This approximation is similar in spirit to the fixed-node approximation in real-space DMC but applied in the Hilbert space where the AFQMC random walk occurs.

Theoretical Background
----------------------

... coming soon ...

Input
-----

The input for an AFQMC calculation is fundamentally different to the
input for other real-space algorithms in QMCPACK. The main source of
input comes from the Hamiltonian matrix elements in an appropriate
single particle basis. This must be evaluated by an external code and
saved in a format that QMCPACK can read. More details about file formats
follow. The input file has six basic xml-blocks: ``AFQMCInfo``,
``Hamiltonian``, ``Wavefunction``, ``WalkerSet``, ``Propagator``, and
``execute``. The first five define input structures required for various
types of calculations. The ``execute`` block represents actual
calculations and takes as input the other blocks. Nonexecution blocks
are parsed first, followed by a second pass where execution blocks are
parsed (and executed) in order. :ref:`Listing 51 <Listing 51>` shows an example of a
minimal input file for an AFQMC calculation.
:numref:`table13` shows a brief description of the most
important parameters in the calculation. All xml sections contain a
“name” argument used to identify the resulting object within QMCPACK.
For example, in the example, multiple Hamiltonian objects with different
names can be defined. The one actually used in the calculation is the
one passed to “execute” as ham.

.. code-block::
  :caption: Sample input file for AFQMC.
  :name: Listing 51

  <?xml version="1.0"?>
  <simulation method="afqmc">
    <project id="Carbon" series="0"/>

    <AFQMCInfo name="info0">
      <parameter name="NMO">32</parameter>
      <parameter name="NAEA">16</parameter>
      <parameter name="NAEB">16</parameter>
    </AFQMCInfo>

    <Hamiltonian name="ham0" info="info0">
      <parameter name="filename">../fcidump.h5</parameter>
    </Hamiltonian>

    <Wavefunction name="wfn0" type="MSD" info="info0">
      <parameter name="filetype">ascii</parameter>
      <parameter name="filename">wfn.dat</parameter>
    </Wavefunction>

    <WalkerSet name="wset0">
      <parameter name="walker_type">closed</parameter>
    </WalkerSet>

    <Propagator name="prop0" info="info0">
    </Propagator>

    <execute wset="wset0" ham="ham0" wfn="wfn0" prop="prop0" info="info0">
      <parameter name="timestep">0.005</parameter>
      <parameter name="blocks">10000</parameter>
      <parameter name="nWalkers">20</parameter>
    </execute>

  </simulation>

The following list includes all input sections for AFQMC calculations, along with a detailed explanation of accepted parameters. Since the code is under active development, the list of parameters and their interpretation might change in the future.

``AFQMCInfo``: Input block that defines basic information about the
calculation. It is passed to all other input blocks to propagate the
basic information: ``<AFQMCInfo name="info0">``

-  **NMO**. Number of molecular orbitals, i.e., number of states in the
   single particle basis.

-  **NAEA**. Number of active electrons-alpha, i.e., number of spin-up
   electrons.

-  **NAEB**. Number of active electrons-beta, i.e., number of spin-down
   electrons.

``Hamiltonian``: Controls the object that reads, stores, and manages the
``hamiltonian``.
``<Hamiltonian name="ham0" type="SparseGeneral" info="info0">``

-  **filename**. Name of file with the ``Hamiltonian``. This is a
   required parameter.

-  **cutoff_1bar**. Cutoff applied to integrals during reading. Any term
   in the Hamiltonian smaller than this value is set to zero. (For
   filetype=“hdf5”, the cutoff is applied only to the 2-electron
   integrals). Default: 1e-8

-  **cutoff_decomposition**. Cutoff used to stop the iterative cycle in
   the generation of the Cholesky decomposition of the 2-electron
   integrals. The generation of Cholesky vectors is stopped when the
   maximum error in the diagonal reaches this value. In case of an
   eigenvalue factorization, this becomes the cutoff applied to the
   eigenvalues. Only eigenvalues above this value are kept. Default:
   1e-6

-  **nblocks**. This parameter controls the distribution of the
   2-electron integrals among processors. In the default behavior
   (nblocks=1), all nodes contain the entire list of integrals. If
   nblocks :math:`>` 1, the of nodes in the calculation will be split in
   nblocks groups. Each node in a given group contains the same subset
   of integrals and subsequently operates on this subset during any
   further operation that requires the hamiltonian. The maximum number
   of groups is NMO. Currently only works for filetype=“hdf5” and the
   file must contain integrals. Not yet implemented for input
   hamiltonians in the form of Cholesky vectors or for ASCII input.
   Coming soon! Default: No distribution

-  **printEig**. If “yes”, prints additional information during the
   Cholesky decomposition. Default: no

-  **fix_2eint**. If this is set to “yes”, orbital pairs that are
   found not to be positive definite are ignored in the generation of
   the Cholesky factorization. This is necessary if the 2-electron
   integrals are not positive definite because of round-off errors in
   their generation. Default: no

``Wavefunction``: controls the object that manages the trial
wavefunctions. This block expects a list of xml-blocks defining actual
trial wavefunctions for various roles.
``<Wavefunction name="wfn0" type="MSD/PHMSD" info="info0">``

-  **filename**. Name of file with wavefunction information.

-  **cutoff**. cutoff applied to the terms in the calculation of the
   local energy. Only terms in the Hamiltonian above this cutoff are
   included in the evaluation of the energy. Default: 1e-6

-  **nnodes**. Defines the parallelization of the local energy
   evaluation and the distribution of the ``Hamiltonian`` matrix (not to
   be confused with the list of 2-electron integrals managed by
   ``Hamiltonian``. These are not the same.) If nnodes :math:`>` 1, the
   nodes in the simulation are split into groups of nnodes, each group
   works collectively in the evaluation of the local energy of their
   walkers. This helps distribute the effort involved in the evaluation
   of the local energy among the nodes in the group, but also
   distributes the memory associated with the wavefunction among the
   nodes in the group. Default: No distribution

-  **ndet**. Number of determinants to read from file. Default: Read all
   determinants.

-  **cutoff**. For sparse hamiltoniants, this defines the cutoff applied
   to the half-rotated 2-electron integrals. Default: 0.0

-  **nbatch**. This turns on(>=1)/off(==0) batched calculation of
   density matrices and overlaps. -1 means all the walkers in the batch.
   Default: 0 (CPU) / -1 (GPU)

-  **nbatch_qr**. This turns on(>=1)/off(==0) batched QR calculation. -1
   means all the walkers in the batch. Default: 0 (CPU) / -1 (GPU)

``WalkerSet``: Controls the object that handles the set of walkers.
``<WalkerSet name="wset0">``

-  **walker_type**. Type of walker set: closed or collinear. Default:
   collinear

-  **pop_control**. Population control algorithm. Options: “simple”:
   Uses a simple branching scheme with a fluctuating population. Walkers
   with weight above max_weight are split into multiple walkers of
   weight reset_weight. Walkers with weight below min_weight are killed
   with probability (weight/min_weight); “pair”: Fixed-population
   branching algorithm, based on QWalk’s branching algorithm. Pairs of
   walkers with weight above/below max_weight/min_weight are combined
   into 2 walkers with weights equal to :math:`(w_1+w_2)/2`. The
   probability of replicating walker w1 (larger weight) occurs with
   probability :math:`w_1/(w_1+w_2)`, otherwise walker w2 (lower weight)
   is replicated; “comb”: Fixed-population branching algorithm based on
   the Comb method. Will be available in the next release. Default:
   “pair”

-  **min_weight**. Weight at which walkers are possibly killed (with
   probability weight/min_weight). Default: 0.05

-  **max_weight**. Weight at which walkers are replicated. Default: 4.0

-  **reset_weight**. Weight to which replicated walkers are reset to.
   Default: 1.0

``Propagator``: Controls the object that manages the propagators.
``<Propagator name="prop0" info="info0">``

-  **cutoff**. Cutoff applied to Cholesky vectors. Elements of the
   Cholesky vectors below this value are set to zero. Only meaningful
   with sparse hamiltonians. Default: 1e-6

-  **substractMF**. If “yes”, apply mean-field subtraction based on the
   ImpSamp trial wavefunction. Must set to “no” to turn it off. Default:
   yes

-  **vbias_bound**. Upper bound applied to the vias potential.
   Components of the vias potential above this value are truncated
   there. The bound is currently applied to
   :math:`\sqrt{\tau} v_{bias}`, so a larger value must be used as
   either the time step or the fluctuations increase (e.g. from running
   a larger system or using a poor trial wavefunction). Default: 3.0

-  **apply_constrain**. If “yes”, apply the phaseless constrain to the
   walker propagation. Currently, setting this to “no” produces unknown
   behavior, since free propagation algorithm has not been tested.
   Default: yes

-  **hybrid**. If “yes”, use hybrid propagation algorithm. This
   propagation scheme doesn’t use the local energy during propagation,
   leading to significant speed ups when its evaluation cost is high.
   The local energy of the ImpSamp trial wavefunction is never
   evaluated. To obtain energy estimates in this case, you must define
   an Estimator xml-block with the ``Wavefunction`` block. The local
   energy of this trial wavefunction is evaluated and printed. It is
   possible to use a previously defined trial wavefunction in the
   Estimator block, just set its “name” argument to the name of a
   previously defined wavefunction. In this case, the same object is
   used for both roles. Default: no

-  **nnodes**. Controls the parallel propagation algorithm. If nnodes
   :math:`>` 1, the nodes in the simulation are split into groups of
   nnodes nodes, each group working collectively to propagate their
   walkers. Default: 1 (Serial algorithm)

-  **nbatch**. This turns on(>=1)/off(==0) batched calculation of
   density matrices and overlaps. -1 means all the walkers in the batch.
   Default: 0 (CPU) / -1 (GPU)

-  **nbatch_qr**. This turns on(>=1)/off(==0) batched QR
   calculation. -1 means all the walkers in the batch. Default: 0 (CPU)
   / -1 (GPU)

``execute``: Defines an execution region.
``<execute wset="wset0" ham="ham0" wfn="wfn0" prop="prop0" info="info0">``

- **nWalkers**. Initial number of walkers per core group (see
  ncores). This sets the number of walkers for a given gorup of
  “ncores" on a node; the total number of walkers in the simulation
  depends on the total number of nodes and on the total number of
  cores on a node in the following way:
  :math:`\#_walkers_total = nWalkers * \#_nodes * \#_cores_total / ncores`.
  Default: 5

- **timestep**. Time step in 1/a.u.
  Default: 0.01

- **blocks**. Number of blocks. Slow operations occur once per block
  (e.g., write to file, slow observables, checkpoints),
  Default: 100

- **step**. Number of steps within a block. Operations that occur at
  the step level include load balance, orthogonalization, branching,
  etc.
  Default: 1

- **substep**. Number of substeps within a step. Only walker
  propagation occurs in a substep.
  Default: 1

- **ortho**. Number of steps between orthogonalization. Default: 1

- **ncores**. Number of nodes in a task group. This number defines the
  number of cores on a node that share the parallel work associated
  with a distributed task. This number is used in the ``Wavefunction``
  and ``Propagator`` task groups. The walker sets are shares by the
  ncores on a given node in the task group.

- **checkpoint**. Number of blocks between checkpoint files are
  generated. If a value smaller than 1 is given, no file is generated.
  If **hdf_write_file** is not set, a default name is used. **Default:
  0**

- **hdf_write_file**. If set (and checkpoint>0), a checkpoint file with
  this name will be written.

- **hdf_read_file**. If set, the simulation will be restarted from
  the given file.

Within the ``Estimators`` xml block has an argument **name**: the type
of estimator we want to measure. Currently available estimators include:
“basic”, “energy”, “mixed_one_rdm”, and “back_propagation”.

The basic estimator has the following optional parameters:

-  **timers**. print timing information. Default: true

The back_propagation estimator has the following parameters:

-  **ortho**. Number of back-propagation steps between
   orthogonalization. Default: 10

-  **nsteps**. Maximum number of back-propagation steps. Default: 10

-  **naverages**. Number of back propagation calculations to perform.
   The number of steps will be chosed equally distributed in the range
   0,nsteps. Default: 1

-  **block_size**. Number of blocks to use in the internal average of
   the back propagated estimator. This is used to block data and reduce
   the size of the output. Default: 1

-  **nskip**. Number of blocks to skip at the start of the calculation
   for equilibration purposes. Default: 0

File formats
------------

QMCPACK offers three factorization approaches which are appropriate in different settings. The most generic approach implemented
is based on the modified-Cholesky
factorization :cite:`BeebeCholesky1977,KochCholesky2003,AquilanteMOLCAS2009,PurwantoCa2011,PurwantoDownfolding2013` of the ERI
tensor:

.. math::
  :label: eq58

  v_{pqrs} = V_{(pr),(sq)} \approx \sum_n^{N_\mathrm{chol}} L_{pr,n} L^{*}_{sq,n},

where the sum is truncated at :math:`N_{\mathrm{chol}} = x_c M`,
:math:`x_c` is typically between :math:`5` and :math:`10`, :math:`M` is
the number of basis functions and we have assumed that the
single-particle orbitals are in general complex. The storage requirement
is thus naively :math:`\mathcal{O}(M^3)`. Note we follow the usual
definition of :math:`v_{pqrs} = \langle pq | rs \rangle = (pr|qs)`. With
this form of factorization QMCPACK allows for the integrals to be stored
in either dense or sparse format.

The dense case is the simplest and is only implemented for Hamiltonians
with *real* integrals (and basis functions, i.e. not the homegeneous
electron gas which has complex orbitals but real integrals). The file
format is given as follows:

.. code-block::
  :caption: Sample Dense Cholesky QMCPACK Hamtiltonian.
  :name: Listing 52

  $ h5dump -n afqmc.h5
  HDF5 "afqmc.h5" {
      FILE_CONTENTS {
          group      /
          group      /Hamiltonian
          group      /Hamiltonian/DenseFactorized
          dataset    /Hamiltonian/DenseFactorized/L
          dataset    /Hamiltonian/dims
          dataset    /Hamiltonian/hcore
          dataset    /Hamiltonian/Energies
      }
  }

where the datasets are given by the following

-  ``/Hamiltonian/DenseFactorized/L`` Contains the :math:`[M^2,N_\mathrm{nchol}]` dimensional matrix
   representatation of :math:`L_{pr,n}`.

-  ``/Hamiltonian/dims`` Descriptor array of length 8 containing
   :math:`[0,0,0,M,N_\alpha,N_\beta,0,N_\mathrm{nchol}]`. Note that
   :math:`N_\alpha` and :math:`N_\beta` are somewhat redundant and will
   be read from the input file and wavefunction. This allows for the
   Hamiltonian to be used with different (potentially spin polarized)
   wavefunctions.

-  ``/Hamiltonian/hcore`` Contains the :math:`[M,M]` dimensional one-body Hamiltonian matrix
   elements :math:`h_{pq}`.

-  ``/Hamiltonian/Energies`` Array containing :math:`[E_{II}, E_{\mathrm{core}}]`.
   :math:`E_{II}` should contain ion-ion repulsion energy and any
   additional constant terms which have to be added to the total energy.
   :math:`E_{\mathrm{core}}` is deprecated and not used.

Typically the Cholesky matrix is sparse, particularly if written in the
non-orthogonal AO basis (not currently supported in QMCPACK). In this
case only a small number of non-zero elements (denoted :math:`nnz`
below) need to be stored which can reduce the memory overhead
considerably. Internally QMCPACK stores this matrix in the CSR format,
and the HDF5 file format is reflective of this. For large systems and,
more generally when running in parallel, it is convenient to chunk the
writing/reading of the Cholesky matrix into blocks of size
:math:`[M^2,\frac{N_{\mathrm{chol}}}{N_{\mathrm{blocks}}}]` (if
interpreted as a dense array). This is achieved by writing these blocks
to different data sets in the file. For the sparse case the Hamtiltonian
file format is given as follows:

.. code-block::
  :caption: Sample Sparse Cholesky QMCPACK Hamtiltonian.
  :name: Listing 53

  $ h5dump -n afqmc.h5
  HDF5 "afqmc.h5" {
      FILE_CONTENTS {
          group      /
          group      /Hamiltonian
          group      /Hamiltonian/Factorized
          dataset    /Hamiltonian/Factorized/block_sizes
          dataset    /Hamiltonian/Factorized/index_0
          dataset    /Hamiltonian/Factorized/vals_0
          dataset    /Hamiltonian/ComplexIntegrals
          dataset    /Hamiltonian/dims
          dataset    /Hamiltonian/hcore
          dataset    /Hamiltonian/Energies
      }
  }

-  ``/Hamiltonian/Factorized/block_sizes`` Contains the number of elements in each block of the sparse
   representation of the Cholesky matrix :math:`L_{pr,n}`. In this case
   there is 1 block.

-  ``/Hamiltonian/Factorized/index_0`` :math:`[2\times nnz]` dimensional array, containing the indices of
   the non-zero values of :math:`L_{ik,n}`. The row indices are stored
   in the even entries, and the column indices in the odd entries.

-  ``/Hamiltonian/Factorized/vals_0`` :math:`[nnz]` length array containing non-zero values of
   :math:`L_{pr,n}` for chunk 0.

-  ``/Hamiltonian/dims`` Descriptor array of length 8 containing
   :math:`[0,nnz,N_{\mathrm{block}},M,N_\alpha,N_\beta,0,N_\mathrm{nchol}]`.

-  ``/Hamiltonian/ComplexIntegrals`` Length 1 array that specifies if integrals are complex valued. 1
   for complex integrals, 0 for real integrals.

-  ``/Hamiltonian/hcore`` Contains the :math:`[M,M]` dimensional one-body Hamiltonian matrix
   elements :math:`h_{pq}`. Due to its small size this is written as a
   dense 2D-array.

-  ``/Hamiltonian/Energies`` Array containing :math:`[E_{II}, E_{\mathrm{core}}]`.
   :math:`E_{II}` should contain ion-ion repulsion energy and any
   additional constant terms which have to be added to the total energy.
   :math:`E_{\mathrm{core}}` is deprecated and not used.

To reduce the memory overhead of storing the three-index tensor we recently adapted the
tensor-hypercontraction :cite:`HohensteinTHCI2012,ParrishTHCII2012,HohensteinTHCIII2012` (THC) approach for use in AFQMC\cite{MaloneISDF2019}. Within the THC approach we
can approximate the orbital products entering the ERIs as

.. math::
  :label: eq59

  \varphi^{*}_p(\mathbf{r})\varphi_r(\mathbf{r}) \approx \sum_\mu^{N_\mu} \zeta_\mu(\mathbf{r}) \varphi^*_p(\mathbf{r}_\mu)\varphi_r(\mathbf{r}_\mu),

where :math:`\varphi_p(\mathbf{r})` are the one-electron orbitals and
:math:`\mathbf{r}_\mu` are a set of specially selected interpolating
points, :math:`\zeta_\mu(\mathbf{r})` are a set of interpolating vectors
and :math:`N_\mu = x_\mu M`. We can then write the ERI tensor as a
product of rank-2 tensors

.. math::
  :label: eq60

  v_{pqrs} \approx \sum_{\mu\nu} \varphi^{*}_p(\mathbf{r}_\mu)\varphi_r(\mathbf{r}_\mu) M_{\mu\nu} \varphi^{*}_q(\mathbf{r}_\nu)\varphi_s(\mathbf{r}_\nu),

where

.. math::
  :label: eq61

  M_{\mu\nu} = \int d\mathbf{r}d\mathbf{r}' \zeta_\mu(\mathbf{r})\frac{1}{|\mathbf{r}-\mathbf{r}'|}\zeta^{*}_\nu(\mathbf{r}').

We also require the half-rotated versions of these quantities which live
on a different set of :math:`\tilde{N}_\mu` interpolating points
:math:`\tilde{\mathbf{r}}_\mu` (see :cite:`MaloneISDF2019`). The file format for THC
factorization is as follows:

.. code-block::
  :caption: Sample Sparse Cholesky QMCPACK Hamtiltonian.
  :name: Listing 54

  $ h5dump -n afqmc.h5
  HDF5 "afqmc.h5" {
      FILE_CONTENTS {
          group      /
          group      /Hamiltonian
          group      /Hamiltonian/THC
          dataset    /Hamiltonian/THC/Luv
          dataset    /Hamiltonian/THC/Orbitals
          dataset    /Hamiltonian/THC/HalfTransformedMuv
          dataset    /Hamiltonian/THC/HalfTransformedFullOrbitals
          dataset    /Hamiltonian/THC/HalfTransformedOccOrbitals
          dataset    /Hamiltonian/THC/dims
          dataset    /Hamiltonian/ComplexIntegrals
          dataset    /Hamiltonian/dims
          dataset    /Hamiltonian/hcore
          dataset    /Hamiltonian/Energies
      }
  }

-  ``/Hamiltonian/THC/Luv`` Cholesky factorization of the :math:`M_{\mu\nu}` matrix given in :eq:`eq61`.

-  ``/Hamiltonian/THC/Orbitals`` :math:`[M,N_\mu]` dimensional array of orbitals evaluated at chosen
   interpolating points :math:`\varphi_i(\mathbf{r}_\mu)`.

-  ``/Hamiltonian/THC/HalfTransformedMuv`` :math:`[\tilde{N}_\mu,\tilde{N}_\mu]` dimensional array containing
   half-transformed :math:`\tilde{M}_{\mu\nu}`.

-  ``/Hamiltonian/THC/HalfTransformedFullOrbitals`` :math:`[M,\tilde{N}_\mu]` dimensional array containing orbital set
   computed at half-transformed interpolating points
   :math:`\varphi_i(\tilde{\mathbf{r}}_\mu)`.

-  ``/Hamiltonian/THC/HalfTransformedOccOrbitals`` :math:`[N_\alpha+N_\beta,\tilde{N}_\mu]` dimensional array
   containing half-rotated orbital set computed at half-transformed
   interpolating points
   :math:`\varphi_a(\tilde{\mathbf{r}}_\mu) = \sum_{p} A_{pa}^* \varphi_{p}(\tilde{\mathbf{r}}_\mu)`,
   where :math:`\mathbf{A}` is the Slater-Matrix of the (currently
   single-determinant) trial wavefunction.

-  ``/Hamiltonian/THC/dims`` Descriptor array containing :math:`[M, N_\mu, \tilde{N}_\mu]`.

-  ``/Hamiltonian/ComplexIntegrals`` Length 1 array that specifies if integrals are complex valued. 1
   for complex integrals, 0 for real integrals.

-  ``/Hamiltonian/dims`` Descriptor array of length 8 containing
   :math:`[0,0,0,M,N_\alpha,N_\beta,0,0]`.

-  ``/Hamiltonian/hcore`` Contains the :math:`[M,M]` dimensional one-body Hamiltonian matrix
   elements :math:`h_{ij}`.

-  ``/Hamiltonian/Energies`` Array containing :math:`[E_{II}, E_{\mathrm{core}}]`.
   :math:`E_{II}` should contain ion-ion repulsion energy and any
   additional constant terms which have to be added to the total energy
   (such as the electron-electron interaction Madelung contribution of
   :math:`\frac{1}{2} N \xi )`. :math:`E_{\mathrm{core}}` is deprecated
   and not used.

Finally, we have implemented an explicitly :math:`k`-point dependent factorization for periodic systems :cite:`MottaKPoint2019,MaloneGPU2020`

.. math::
  :label: eq62

  v_{pqrs} = \sum_{\substack{n\textbf{Q}\textbf{k}\textbf{k}' \\ pqrs\sigma\sigma'}} L^{\textbf{Q},\textbf{k}}_{pr,n} {L^{\textbf{Q},\textbf{k}'}_{sq,n}}^{*}

where :math:`\textbf{k}`, :math:`\textbf{k}'` and :math:`\textbf{Q}` are
vectors in the first Brillouin zone. The one-body Hamiltonian is block
diagonal in :math:`\textbf{k}` and in :eq:`eq62` we have used
momentum conservation
:math:`(\textbf{k}_p - \textbf{k}_r + \textbf{k}_q - \textbf{k}_s) = \textbf{G}`
with :math:`\textbf{G}` being some vector in the reciprocal lattice of
the simulation cell. The convention for the Cholesky matrix
:math:`L^{\textbf{Q},\textbf{k}}_{pr,\gamma}` is as follows:
:math:`\textbf{k}_r = \textbf{k}_p - \textbf{Q}`, so the vector
:math:`\textbf{k}` labels the *k*-point of the first band index,
:math:`\textit{p}`, while the *k*-point vector of the second band index,
:math:`\textit{r}`, is given by :math:`\textbf{k} - \textbf{Q}`.
Electron repulsion integrals at different :math:`\textbf{Q}` vectors are
zero by symmetry, resulting in a reduction in the number of required
:math:`\mathbf{Q}` vectors. For certain :math:`\textbf{Q}` vectors that
satisfy :math:`\textbf{Q} \ne -\textbf{Q}` (this is not satisfied at the
origin and at high symmetry points on the edge of the 1BZ), we have
:math:`{L^{\textbf{Q},\textbf{k}}_{sq,\gamma}}^{*} = {L^{-\textbf{Q},\textbf{k}-\textbf{Q}}_{qs,\gamma}}`,
which requires us to store Cholesky vectors for either one of the
:math:`(\textbf{Q},-\textbf{Q})` pair, but not both.

In what follows let :math:`m_{\mathbf{k}}` denote the number of basis
functions for basis functions of a given :math:`k`-point (these can in
principle differ for different :math:`k`-points due to linear
dependencies), :math:`n^{\alpha}_{\mathbf{k}}` the number of
:math:`\alpha` electrons in a given :math:`k`-point and
:math:`n_{\mathrm{chol}}^{\mathbf{Q}_n}` the number of Cholesky vectors
for momentum transfer :math:`\mathbf{Q}_n`. The file format for this
factorization is as follows (for a :math:`2\times2\times2`
:math:`k`-point mesh, for denser meshes generally there will be far
fewer symmetry inequivalent momentum transfer vectors than there are
:math:`k`-points):

.. code-block::
  :caption: Sample Dense :math:`k`-point dependent Cholesky QMCPACK Hamtiltonian.
  :name: Listing 55

  $ h5dump -n afqmc.h5
  HDF5 "afqmc.h5" {
      FILE_CONTENTS {
          group      /
          group      /Hamiltonian
          group      /Hamiltonian/KPFactorized
          dataset    /Hamiltonian/KPFactorized/L0
          dataset    /Hamiltonian/KPFactorized/L1
          dataset    /Hamiltonian/KPFactorized/L2
          dataset    /Hamiltonian/KPFactorized/L3
          dataset    /Hamiltonian/KPFactorized/L4
          dataset    /Hamiltonian/KPFactorized/L5
          dataset    /Hamiltonian/KPFactorized/L6
          dataset    /Hamiltonian/KPFactorized/L7
          dataset    /Hamiltonian/NCholPerKP
          dataset    /Hamiltonian/MinusK
          dataset    /Hamiltonian/NMOPerKP
          dataset    /Hamiltonian/QKTok2
          dataset    /Hamiltonian/H1_kp0
          dataset    /Hamiltonian/H1_kp1
          dataset    /Hamiltonian/H1_kp2
          dataset    /Hamiltonian/H1_kp3
          dataset    /Hamiltonian/H1_kp4
          dataset    /Hamiltonian/H1_kp5
          dataset    /Hamiltonian/H1_kp6
          dataset    /Hamiltonian/H1_kp7
          dataset    /Hamiltonian/ComplexIntegrals
          dataset    /Hamiltonian/KPoints
          dataset    /Hamiltonian/dims
          dataset    /Hamiltonian/Energies
      }
  }

-  ``/Hamiltonian/KPFactorized/L[n]`` This series of datasets store elements of the Cholesky tensors
   :math:`L[\mathbf{Q}_n,\mathbf{k},pr,n]`. Each data set is of
   dimension
   :math:`[N_k,m_{\mathbf{k}}\times m_{\mathbf{k}'},n^{\mathbf{Q}_n}_\mathrm{chol}]`,
   where, again, :math:`k` is the :math:`k`-point associated with basis
   function :math:`p`, the :math:`k`-point of basis function :math:`r`
   is defined via the mapping ``QKtok2``.

-  ``/Hamiltonian/NCholPerKP`` :math:`N_k` length array giving number of Cholesky vectors per
   :math:`k`-point.

-  ``/Hamiltonian/MinusK``: :math:`N_k` length array mapping a
   :math:`k`-point to its inverse: :math:`\mathbf{k}_i+`\ ``MinusK[i]``
   :math:`= \mathbf{0} \mod \mathbf{G}`.

-  ``/Hamiltonian/NMOPerKP``: :math:`N_k` length array listing number of
   basis functions per :math:`k`-point.

-  ``/Hamiltonian/QKTok2``: :math:`[N_k,N_k]` dimensional array.
   ``QKtok2[i,j]`` yields the :math:`k` point index satisfying
   :math:`\mathbf{k}=\mathbf{Q}_i-\mathbf{k}_j+\mathbf{G}`.

-  ``/Hamiltonian/dims``: Descriptor array of length 8 containing
   :math:`[0,0,0,M,N_\alpha,N_\beta,0,0]`.

-  ``/Hamiltonian/H1_kp[n]`` Contains the :math:`[m_{\mathbf{k}_n},m_{\mathbf{k}_n}]`
   dimensional one-body Hamiltonian matrix elements
   :math:`h_{(\mathbf{k}_{n}p)(\mathbf{k}_{n}q)}`.

-  ``/Hamiltonian/ComplexIntegrals`` Length 1 array that specifies if integrals are complex valued. 1
   for complex integrals, 0 for real integrals.

-  ``/Hamiltonian/KPoints`` :math:`[N_k,3]` Dimensional array containing :math:`k`-points used to
   sample Brillouin zone.

-  ``/Hamiltonian/dims`` Descriptor array of length 8 containing
   :math:`[0,0,N_k,M,N_\alpha,N_\beta,0,N_\mathrm{nchol}]`. Note that
   :math:`M` is the total number of basis functions, i.e.
   :math:`M=\sum_\mathbf{k} m_\mathbf{k}`, and likewise for the number
   of electrons.

-  ``/Hamiltonian/Energies`` Array containing :math:`[E_{II}, E_{\mathrm{core}}]`.
   :math:`E_{II}` should contain ion-ion repulsion energy and any
   additional constant terms which have to be added to the total energy
   (such as the electron-electron interaction Madelung contribution of
   :math:`\frac{1}{2} N \xi )`. :math:`E_{\mathrm{core}}` is deprecated
   and not used.

Complex integrals should be written as an array with an additional dimension, e.g., a 1D array should be written as a 2D array with ``array_hdf5[:,0]=real(1d_array)`` and ``array_hdf5[:,1]=imag(1d_array)``. The functions ``afqmctools.utils.misc.from_qmcpack_complex`` and ``afqmctools.utils.misc.to_qmcpack_complex`` can be used to transform qmcpack format to complex valued numpy arrays of the appropriate shape and vice versa.

Finally, if using external tools to generate this file format, we provide a sanity checker script in ``utils/afqmctools/bin/test_afqmc_input.py`` which will raise errors if the format does not conform to what is being used internally.

Advice/Useful Information
-------------------------

AFQMC calculations are computationally expensive and require some care to obtain reasonable performance.
The following is a growing list of useful advice for new users, followed by a sample input for a large calculation.

-  Generate Cholesky-decomposed integrals with external codes instead of
   the 2-electron integrals directly. The generation of the Cholesky
   factorization is faster and consumes less memory.

-  Use the hybrid algorithm for walker propagation. Set steps/substeps
   to adequate values to reduce the number of energy evaluations. This
   is essential when using large multideterminant expansions.

-  Adjust cutoffs in the wavefunction and propagator bloxks until
   desired accuracy is reached. The cost of the calculation will depend
   on these cutoffs.

-  Adjust ncores/nWalkers to obtain better efficiency. Larger nWalkers
   will lead to more efficient linear algebra operations but will
   increase the time per step. Larger ncores will reduce the time per
   step but will reduce efficiency because of inefficiencies in the
   parallel implementation. For large calculations, values between 6–12
   for both quantities should be reasonable, depending on architecture.

.. code-block::
  :caption: Example of sections of an AFQMC input file for a large calculation.
  :name: Listing 56

  ...

    <Hamiltonian name="ham0" type="SparseGeneral" info="info0">
      <parameter name="filename">fcidump.h5</parameter>
      <parameter name="cutoff_1bar">1e-6</parameter>
      <parameter name="cutoff_decomposition">1e-5</parameter>
    </Hamiltonian>

    <Wavefunction name="wfn0" type="MSD" info="info0">
      <parameter name="filetype">ascii</parameter>
      <parameter name="filename">wfn.dat</parameter>
    </Wavefunction>

    <WalkerSet name="wset0">
      <parameter name="walker_type">closed</parameter>
    </WalkerSet>

    <Propagator name="prop0" info="info0">
      <parameter name="hybrid">yes</parameter>
    </Propagator>

    <execute wset="wset0" ham="ham0" wfn="wfn0" prop="prop0" info="info0">
      <parameter name="ncores">8</parameter>
      <parameter name="timestep">0.01</parameter>
      <parameter name="blocks">10000</parameter>
      <parameter name="steps">10</parameter>
      <parameter name="substeps">5</parameter>
      <parameter name="nWalkers">8</parameter>
      <parameter name="ortho">5</parameter>
    </execute>

.. centered:: ``afqmc`` method

parameters in ``AFQMCInfo``

.. _table13:
.. table::

  +----------+--------------+---------------------------+-------------+-----------------------------------------+
  | **Name** | **Datatype** | **Values**                | **Default** | **Description**                         |
  +==========+==============+===========================+=============+=========================================+
  | ``NMO``  | integer      | :math:`\geq  0`           | no          | Number of molecular orbitals            |
  +----------+--------------+---------------------------+-------------+-----------------------------------------+
  | ``NAEA`` | integer      | :math:`\geq  0`           | no          | Number of active electrons of spin-up   |
  +----------+--------------+---------------------------+-------------+-----------------------------------------+
  | ``NAEB`` | integer      | :math:`\geq  0`           | no          | Number of active electrons of spin-down |
  +----------+--------------+---------------------------+-------------+-----------------------------------------+

parameters in ``Hamiltonian``

+--------------+--------------+------------+-------------+-------------------------------------+
| **Name**     | **Datatype** | **Values** | **Default** | **Description**                     |
+==============+==============+============+=============+=====================================+
| ``info``     | argument     |            |             | Name of ``AFQMCInfo`` block         |
+--------------+--------------+------------+-------------+-------------------------------------+
| ``filename`` | string       |            | no          | Name of file with the hamiltonian   |
+--------------+--------------+------------+-------------+-------------------------------------+
| ``filetype`` | string       | hdf5       | yes         | Native HDF5-based format of QMCPACK |
+--------------+--------------+------------+-------------+-------------------------------------+

parameters in ``Wavefunction``

+--------------+--------------+-------------+-------------+--------------------------------------------------------------------+
| **Name**     | **Datatype** | **Values**  | **Default** | **Description**                                                    |
+==============+==============+=============+=============+====================================================================+
| ``info``     | argument     |             |             | name of ``AFQMCInfo`` block                                        |
+--------------+--------------+-------------+-------------+--------------------------------------------------------------------+
| ``type``     | argument     | MSD, PHMSD  | no          | Linear combination of (assumed non-orthogonal) Slater determinants |
+--------------+--------------+-------------+-------------+--------------------------------------------------------------------+
| ``filetype`` | string       | ascii, hdf5 | no          | CI-type multi-determinant wave function                            |
+--------------+--------------+-------------+-------------+--------------------------------------------------------------------+

parameters in ``WalkerSet``

+------------------+--------------+------------+-------------+--------------------------------------------------------+
| **Name**         | **Datatype** | **Values** | **Default** | **Description**                                        |
+==================+==============+============+=============+========================================================+
| ``walker_type``  | string       | collinear  | yes         | Request a collinear walker set.                        |
+------------------+--------------+------------+-------------+--------------------------------------------------------+
|                  |              | closed     | no          | Request a closed shell (doubly-occupied) walker set.   |
+------------------+--------------+------------+-------------+--------------------------------------------------------+

parameters in ``Propagator``

+------------+--------------+------------+-------------+------------------------------------------------+
| **Name**   | **Datatype** | **Values** | **Default** | **Description**                                |
+============+==============+============+=============+================================================+
| ``type``   | argument     | afqmc      | afqmc       | Type of propagator                             |
+------------+--------------+------------+-------------+------------------------------------------------+
| ``info``   | argument     |            |             | Name of ``AFQMCInfo`` block                    |
+------------+--------------+------------+-------------+------------------------------------------------+
| ``hybrid`` | string       | yes        |             | Use hybrid propagation algorithm.              |
+------------+--------------+------------+-------------+------------------------------------------------+
|            |              | no         |             | Use local energy based propagation algorithm.  |
+------------+--------------+------------+-------------+------------------------------------------------+

parameters in ``execute``

+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| **Name**     | **Datatype** | **Values**              | **Default** | **Description**                                   |
+==============+==============+=========================+=============+===================================================+
| ``wset``     | argument     |                         |             |                                                   |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``ham``      | argument     |                         |             |                                                   |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``wfn``      | argument     |                         |             |                                                   |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``prop``     | argument     |                         |             |                                                   |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``info``     | argument     |                         |             | Name of ``AFQMCInfo`` block                       |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``nWalkers`` | integer      | :math:`\geq 0`          | 5           | Initial number of walkers per task group          |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``timestep`` | real         | :math:`> 0`             | 0.01        | Time step in 1/a.u.                               |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``blocks``   | integer      | :math:`\geq 0`          | 100         | Number of blocks                                  |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``step``     | integer      | :math:`> 0`             | 1           | Number of steps within a block                    |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``substep``  | integer      | :math:`> 0`             | 1           | Number of substeps within a step                  |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+
| ``ortho``    | integer      | :math:`> 0`             | 1           | Number of steps between walker orthogonalization. |
+--------------+--------------+-------------------------+-------------+---------------------------------------------------+


.. _pyscf:

Using PySCF to generate integrals and trial wavefunctions for AFQMC
------------------------------------------------------------------

.. bibliography:: /bibs/afqmc.bib
