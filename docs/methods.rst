.. _qmcmethods:

Quantum Monte Carlo Methods
===========================

``qmc`` factory element:

  +-----------------+----------------------+
  | Parent elements | ``simulation, loop`` |
  +-----------------+----------------------+
  | type selector   | ``method`` attribute |
  +-----------------+----------------------+

  type options:

  +--------+-----------------------------------------------+
  | vmc    | Variational Monte Carlo                       |
  +--------+-----------------------------------------------+
  | linear | Wavefunction optimization with linear method  |
  +--------+-----------------------------------------------+
  | dmc    | Diffusion Monte Carlo                         |
  +--------+-----------------------------------------------+
  | rmc    | Reptation Monte Carlo                         |
  +--------+-----------------------------------------------+

  shared attributes:

  +----------------+--------------+--------------+-------------+---------------------------------+
  | **Name**       | **Datatype** | **Values**   | **Default** | **Description**                 |
  +================+==============+==============+=============+=================================+
  | ``method``     | text         | listed above | invalid     | QMC driver                      |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``move``       | text         | pbyp, alle   | pbyp        | Method used to move electrons   |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``gpu``        | text         | yes/no       | dep.        | Use the GPU                     |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``trace``      | text         |              | no          | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``profiling``  | text         | yes/no       | no          | Activate resume/pause control   |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``checkpoint`` | integer      | -1, 0, n     | -1          | Checkpoint frequency            |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``record``     | integer      | n            | 0           | Save configuration ever n steps |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``target``     | text         |              |             | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``completed``  | text         |              |             | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+
  | ``append``     | text         | yes/no       | no          | ???                             |
  +----------------+--------------+--------------+-------------+---------------------------------+

Additional information:

-  ``move``: There are two ways to move electrons. The more used method
   is the particle-by-particle move. In this method, only one electron
   is moved for acceptance or rejection. The other method is the
   all-electron move; namely, all the electrons are moved once for
   testing acceptance or rejection.

-  ``gpu``: When the executable is compiled with CUDA, the target
   computing device can be chosen by this switch. With a regular
   CPU-only compilation, this option is not effective.

-  ``profiling``: Performance profiling tools by default profile complete application executions.
   This is largely unnecessary if the focus is a QMC section instead of any initialization
   and additional QMC sections for equilibrating walkers.
   Setting this flag to ``yes`` for the QMC sections of interest and starting the tool with
   data collection paused from the beginning help reducing the profiling workflow
   and amount of collected data. Additional restriction may be imposed by profiling tools.
   For example, NVIDIA profilers can only be turned on and off once and thus only the first QMC
   section with ``profiling="yes"`` will be profiled.
   VTune instead allows pause and resume for unlimited times and thus multiple selected QMC sections
   can be profiled in a single run.

-  ``checkpoint``: This enables and disables checkpointing and
   specifying the frequency of output. Possible values are:

   - **[-1]** No checkpoint (default setting).

   - **[0]** Write the checkpoint files after the completion of the QMC section.

   - **[n]** Write the checkpoint files after every :math:`n` blocks, and also at the end of the QMC section.

The particle configurations are written to a ``.config.h5`` file.

.. code-block::
  :caption: The following is an example of running a simulation that can be restarted.
  :name: Listing 42

  <qmc method="dmc" move="pbyp"  checkpoint="0">
    <parameter name="timestep">         0.004  </parameter>
    <parameter name="blocks">           100   </parameter>
    <parameter name="steps">            400    </parameter>
  </qmc>

The checkpoint flag instructs QMCPACK to output walker configurations.  This also
works in VMC.  This outputs an h5 file with the name ``projectid.run-number.config.h5``.
Check that this file exists before attempting a restart.

To continue a run, specify the ``mcwalkerset`` element before your VMC/DMC block:

.. code-block::
  :caption: Restart (read walkers from previous run).
  :name: Listing 43

  <mcwalkerset fileroot="BH.s002" version="0 6" collected="yes"/>
   <qmc method="dmc" move="pbyp"  checkpoint="0">
     <parameter name="timestep">         0.004  </parameter>
     <parameter name="blocks">           100   </parameter>
     <parameter name="steps">            400    </parameter>
   </qmc>

``BH`` is the project id, and ``s002`` is the calculation number to read in the walkers from the previous run.

In the project id section, make sure that the series number is different from any existing ones to avoid overwriting them.


.. _batched_drivers:

Batched drivers
---------------

The batched drivers introduce a new concept, "crowd", as a sub-organization of walker population.
A crowd is a subset of the walkers that are operated on as as single batch.
Walkers within a crowd operate their computation in lock-step, which helps the GPU efficiency.
Walkers in different crowds remain fully asynchronous unless operations involving the full population are needed.
With this flexible batching capability the new drivers are capable of delivering maximal performance on given hardware.
In the new driver design, all the batched API calls may fallback to an existing single walker implementation.
Consequently, batched drivers allow mixing and matching CPU-only and GPU-accelerated features
in a way that is not feasible with the legacy GPU implementation.



.. _transition_guide:

Transition from classic drivers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Available drivers in batched versions are ``vmc``, ``dmc`` and ``linear``.
There are notable changes in the driver input section when moving from classic drivers to batched drivers:

  - ``walkers`` is not supported in any batched driver inputs.
    Instead, ``walkers_per_rank`` and ``total_walkers`` specify the population at the start of a driver run.

  - ``crowds`` can added in batched drivers to specify the number of crowds.

  - If a classic driver input section contains ``walkers`` equals 1, the same effect can be achieved by
    omitting the specification of ``walkers_per_rank``, ``total_walkers`` or ``crowds`` in batched drivers.

  - The ``walkers_per_rank``, ``total_walkers`` or ``crowds`` parameters are optional.
    See driver-specific parameter additional information below about default values.

  - When running on GPUs, tuning ``walkers_per_rank`` or ``total_walkers`` is likely needed to maximize GPU throughput,
    just like tuning ``walkers`` in the classic drivers.

  - Only particle-by-particle move is supported. No all-particle move support.

  - During development the new drivers had separate names (``vmc_batch``, ``dmc_batch``, and ``linear_batch``).  The use of separate names has been replaced by the ``driver_version`` parameter in the ``project`` section.

.. _vmc:

Variational Monte Carlo
-----------------------

``vmc`` driver
~~~~~~~~~~~~~~

  parameters:

  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | **Name**                       | **Datatype** | **Values**              | **Default** | **Description**                               |
  +================================+==============+=========================+=============+===============================================+
  | ``walkers``                    | integer      | :math:`> 0`             | dep.        | Number of walkers per MPI task                |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``blocks``                     | integer      | :math:`\geq 0`          | 1           | Number of blocks                              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``steps``                      | integer      | :math:`\geq 0`          | 1           | Number of steps per block                     |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``warmupsteps``                | integer      | :math:`\geq 0`          | 0           | Number of steps for warming up                |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``substeps``                   | integer      | :math:`\geq 0`          | 1           | Number of substeps per step                   |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``usedrift``                   | text         | yes,no                  | yes         | Use the algorithm with drift                  |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``timestep``                   | real         | :math:`> 0`             | 0.1         | Time step for each electron move              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``samples``                    | integer      | :math:`\geq 0`          | 0           | Number of walker samples for DMC/optimization |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``stepsbetweensamples``        | integer      | :math:`> 0`             | 1           | Period of sample accumulation                 |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``samplesperthread``           | integer      | :math:`\geq 0`          | 0           | Number of samples per thread                  |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``storeconfigs``               | integer      | all values              | 0           | Write configurations to files                 |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``blocks_between_recompute``   | integer      | :math:`\geq 0`          | dep.        | Wavefunction recompute frequency              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``spinMass``                   | real         | :math:`> 0`             | 1.0         | Effective mass for spin sampling              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``debug_checks``               | text         | see additional info     | dep.        | Turn on/off additional recompute and checks   |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+

Additional information:

- ``walkers`` The number of walkers per MPI task. The initial default number of \ixml{walkers} is one per OpenMP thread or per MPI
  task if threading is disabled. The number is rounded down to a multiple of the number of threads with a minimum of one per
  thread to ensure perfect load balancing. One walker per thread is created in the event fewer ``walkers`` than threads are
  requested.

- ``blocks`` This parameter is universal for all the QMC
  methods. The MC processes are divided into a number of
  ``blocks``, each containing a number of steps. At the end of each block,
  the statistics accumulated in the block are dumped into files,
  e.g., ``scalar.dat``. Typically, each block should have a sufficient number of steps that the I/O at the end of each block is negligible
  compared with the computational cost. Each block should not take so
  long that monitoring its progress is difficult. There should be a
  sufficient number of ``blocks`` to perform statistical analysis.

- ``warmupsteps`` - ``warmupsteps`` are used only for
  equilibration. Property measurements are not performed during
  warm-up steps.

- ``steps`` - ``steps`` are the number of energy and other property measurements to perform per block.

- ``substeps``  For each substep, an attempt is made to move each of the electrons once only by either particle-by-particle or an
  all-electron move.  Because the local energy is evaluated only at
  each full step and not each substep, ``substeps`` are computationally cheaper
  and can be used to reduce the correlation between property measurements
  at a lower cost.

- ``usedrift`` The VMC is implemented in two algorithms with
  or without drift. In the no-drift algorithm, the move of each
  electron is proposed with a Gaussian distribution. The standard
  deviation is chosen as the time step input. In the drift algorithm,
  electrons are moved by Langevin dynamics.

- ``timestep`` The meaning of time step depends on whether or not
  the drift is used. In general, larger time steps reduce the
  time correlation but might also reduce the acceptance ratio,
  reducing overall statistical efficiency. For VMC, typically the
  acceptance ratio should be close to 50% for an efficient
  simulation.

- ``samples`` Separate from conventional energy and other
  property measurements, samples refers to storing whole electron
  configurations in memory ("walker samples") as would be needed by subsequent
  wavefunction optimization or DMC steps. *A standard VMC run to
  measure the energy does not need samples to be set.*

  .. math::

     \texttt{samples}=
     \frac{\texttt{blocks}\cdot\texttt{steps}\cdot\texttt{walkers}}{\texttt{stepsbetweensamples}}\cdot\texttt{number of MPI tasks}

- ``samplesperthread`` This is an alternative way to set the target amount of samples and can be useful when preparing a stored
  population for a subsequent DMC calculation.

  .. math::

     \texttt{samplesperthread}=
     \frac{\texttt{blocks}\cdot\texttt{steps}}{\texttt{stepsbetweensamples}}

- ``stepsbetweensamples`` Because samples generated by consecutive steps are correlated, having ``stepsbetweensamples`` larger
  than 1 can be used to reduces that correlation. In practice, using larger substeps is cheaper than using ``stepsbetweensamples``
  to decorrelate samples.
  
- ``storeconfigs`` If ``storeconfigs`` is set to a nonzero value, then electron configurations during the VMC run are saved to
  files.

- ``blocks_between_recompute`` Recompute the accuracy critical determinant part of the wavefunction
  from scratch: =1 by default when using mixed precision. =0 (no
  recompute) by default when not using mixed precision. Recomputing
  introduces a performance penalty dependent on system size.

- ``spinMass`` Optional parameter to allow the user to change the rate of spin sampling. If spin sampling is on using ``spinor`` == yes in the electron ParticleSet input,  the spin mass determines the rate
  of spin sampling, resulting in an effective spin timestep :math:`\tau_s = \frac{\tau}{\mu_s}`. The algorithm is described in detail in :cite:`Melton2016-1` and :cite:`Melton2016-2`.

- ``debug_checks`` valid values are 'no', 'all', 'checkGL_after_moves'. If the build type is `debug`, the default value is 'all'. Otherwise, the default value is 'no'.

An example VMC section for a simple VMC run:

::

  <qmc method="vmc" move="pbyp">
    <estimator name="LocalEnergy" hdf5="no"/>
    <parameter name="walkers">    256 </parameter>
    <parameter name="warmupSteps">  100 </parameter>
    <parameter name="substeps">  5 </parameter>
    <parameter name="blocks">  20 </parameter>
    <parameter name="steps">  100 </parameter>
    <parameter name="timestep">  1.0 </parameter>
    <parameter name="usedrift">   yes </parameter>
  </qmc>

Here we set 256 ``walkers`` per MPI, have a brief initial equilibration of 100 ``steps``, and then have 20 ``blocks`` of 100 ``steps`` with 5 ``substeps`` each.

The following is an example of VMC section storing configurations (walker samples) for optimization.

::

  <qmc method="vmc" move="pbyp" gpu="yes">
     <estimator name="LocalEnergy" hdf5="no"/>
     <parameter name="walkers">    256 </parameter>
     <parameter name="samples">    2867200 </parameter>
     <parameter name="stepsbetweensamples">    1 </parameter>
     <parameter name="substeps">  5 </parameter>
     <parameter name="warmupSteps">  5 </parameter>
     <parameter name="blocks">  70 </parameter>
     <parameter name="timestep">  1.0 </parameter>
     <parameter name="usedrift">   no </parameter>
   </qmc>

.. _vmc_batch:

Batched ``vmc`` driver (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  parameters:

  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | **Name**                       | **Datatype** | **Values**              | **Default** | **Description**                                 |
  +================================+==============+=========================+=============+=================================================+
  | ``total_walkers``              | integer      | :math:`> 0`             | 1           | Total number of walkers over all MPI ranks      |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``walkers_per_rank``           | integer      | :math:`> 0`             | 1           | Number of walkers per MPI rank                  |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``crowds``                     | integer      | :math:`> 0`             | dep.        | Number of desynchronized dwalker crowds         |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``blocks``                     | integer      | :math:`\geq 0`          | 1           | Number of blocks                                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``steps``                      | integer      | :math:`\geq 0`          | 1           | Number of steps per block                       |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``warmupsteps``                | integer      | :math:`\geq 0`          | 0           | Number of steps for warming up                  |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``substeps``                   | integer      | :math:`\geq 0`          | 1           | Number of substeps per step                     |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``usedrift``                   | text         | yes,no                  | yes         | Use the algorithm with drift                    |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``timestep``                   | real         | :math:`> 0`             | 0.1         | Time step for each electron move                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``samples`` (not ready)        | integer      | :math:`\geq 0`          | 0           | Number of walker samples for in this VMC run    |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``storeconfigs`` (not ready)   | integer      | all values              | 0           | Write configurations to files                   |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``blocks_between_recompute``   | integer      | :math:`\geq 0`          | dep.        | Wavefunction recompute frequency                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``crowd_serialize_walkers``    | integer      | yes, no                 | no          | Force use of single walker APIs (for testing)   |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``debug_checks``               | text         | see additional info     | dep.        | Turn on/off additional recompute and checks     |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``spin_mass``                  | real         | :math:`\geq 0`          | 1.0         | Effective mass for spin sampling                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``measure_imbalance``          | text         | yes,no                  | no          | Measure load imbalance at the end of each block |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+


Additional information:

- ``crowds`` The number of crowds that the walkers are subdivided into on each MPI rank. If not provided, it is set equal to the number of OpenMP threads.

- ``walkers_per_rank`` The number of walkers per MPI rank. The exact number of walkers will be generated before performing random walking.
  It is not required to be a multiple of the number of OpenMP threads. However, to avoid any idle resources, it is recommended to be at
  least the number of OpenMP threads for pure CPU runs. For GPU runs, a scan of this parameter is necessary to reach reasonable single rank
  efficiency and also get a balanced time to solution.
  If neither ``total_walkers`` nor ``walkers_per_rank`` is provided, ``walkers_per_rank`` is set equal to ``crowds``.

- ``total_walkers`` Total number of walkers over all MPI ranks. if not provided, it is computed as ``walkers_per_rank`` times the number of MPI ranks. If both ``total_walkers`` and ``walkers_per_rank`` are provided, ``total_walkers`` must be equal to ``walkers_per_rank`` times the number MPI ranks.

- ``blocks`` This parameter is universal for all the QMC methods. The MC processes are divided into a number of
  ``blocks``, each containing a number of steps. At the end of each block, the statistics accumulated in the block are dumped into files,
  e.g., ``scalar.dat``. Typically, each block should have a sufficient number of steps that the I/O at the end of each block is negligible
  compared with the computational cost. Each block should not take so long that monitoring its progress is difficult. There should be a
  sufficient number of ``blocks`` to perform statistical analysis.

- ``warmupsteps`` - ``warmupsteps`` are used only for
  equilibration. Property measurements are not performed during
  warm-up steps.

- ``steps`` - ``steps`` are the number of energy and other property measurements to perform per block.

- ``substeps``  For each substep, an attempt is made to move each of the electrons once only by either particle-by-particle or an
  all-electron move.  Because the local energy is evaluated only at
  each full step and not each substep, ``substeps`` are computationally cheaper
  and can be used to de-correlation at a low computational cost.

- ``usedrift`` The VMC is implemented in two algorithms with
  or without drift. In the no-drift algorithm, the move of each
  electron is proposed with a Gaussian distribution. The standard
  deviation is chosen as the time step input. In the drift algorithm,
  electrons are moved by Langevin dynamics.

- ``timestep`` The meaning of time step depends on whether or not
  the drift is used. In general, larger time steps reduce the
  time correlation but might also reduce the acceptance ratio,
  reducing overall statistical efficiency. For VMC, typically the
  acceptance ratio should be close to 50% for an efficient
  simulation.

- ``samples`` (not ready)

- ``storeconfigs`` If ``storeconfigs`` is set to a nonzero value, then electron configurations during the VMC run are saved to
  files.

- ``blocks_between_recompute`` Recompute the accuracy critical determinant part of the wavefunction
  from scratch: =1 by default when using mixed precision. =0 (no
  recompute) by default when not using mixed precision. Recomputing
  introduces a performance penalty dependent on system size.

- ``debug_checks`` valid values are 'no', 'all', 'checkGL_after_load', 'checkGL_after_moves', 'checkGL_after_tmove'. If the build type is `debug`, the default value is 'all'. Otherwise, the default value is 'no'.

- ``spin_mass`` Optional parameter to allow the user to change the rate of spin sampling. If spin sampling is on using ``spinor`` == yes in the electron ParticleSet input,  the spin mass determines the rate
  of spin sampling, resulting in an effective spin timestep :math:`\tau_s = \frac{\tau}{\mu_s}`. The algorithm is described in detail in :cite:`Melton2016-1` and :cite:`Melton2016-2`.

An example VMC section for a simple batched ``vmc`` run:

::

  <qmc method="vmc" move="pbyp">
    <estimator name="LocalEnergy" hdf5="no"/>
    <parameter name="walkers_per_rank">    256 </parameter>
    <parameter name="warmupSteps">  100 </parameter>
    <parameter name="substeps">  5 </parameter>
    <parameter name="blocks">  20 </parameter>
    <parameter name="steps">  100 </parameter>
    <parameter name="timestep">  1.0 </parameter>
    <parameter name="usedrift">   yes </parameter>
  </qmc>

Here we set 256 walkers per MPI rank, have a brief initial equilibration of 100 ``steps``, and then have 20 ``blocks`` of 100 ``steps`` with 5 ``substeps`` each.

.. _optimization:

Wavefunction optimization
-------------------------

Optimizing wavefunction is critical in all kinds of real-space QMC calculations
because it significantly improves both the accuracy and efficiency of computation.
However, it is very difficult to directly adopt deterministic minimization approaches because of the stochastic nature of evaluating quantities with MC.
Thanks to the algorithmic breakthrough during the first decade of this century and the tremendous computer power available,
it is now feasible to optimize tens of thousands of parameters in a wavefunction for a solid or molecule.
QMCPACK has multiple optimizers implemented based on the state-of-the-art linear method.
We are continually improving our optimizers for robustness and friendliness and are trying to provide a single solution.
Because of the large variation of wavefunction types carrying distinct characteristics, using several optimizers might be needed in some cases.
We strongly suggested reading recommendations from the experts who maintain these optimizers.

A typical optimization block looks like the following. It starts with method="linear" and contains three blocks of parameters.

::

  <loop max="10">
   <qmc method="linear" move="pbyp" gpu="yes">
     <!-- Specify the VMC options -->
     <parameter name="walkers">              256 </parameter>
     <parameter name="samples">          2867200 </parameter>
     <parameter name="stepsbetweensamples">    1 </parameter>
     <parameter name="substeps">               5 </parameter>
     <parameter name="warmupSteps">            5 </parameter>
     <parameter name="blocks">                70 </parameter>
     <parameter name="timestep">             1.0 </parameter>
     <parameter name="usedrift">              no </parameter>
     <estimator name="LocalEnergy" hdf5="no"/>
     ...
     <!-- Specify the correlated sampling options and define the cost function -->
     <parameter name="minwalkers">            0.3 </parameter>
          <cost name="energy">               0.95 </cost>
          <cost name="unreweightedvariance"> 0.00 </cost>
          <cost name="reweightedvariance">   0.05 </cost>
     ...
     <!-- Specify the optimizer options -->
     <parameter name="MinMethod">    OneShiftOnly </parameter>
     ...
   </qmc>
  </loop>

  -  Loop is helpful to repeatedly execute identical optimization blocks.

  -  The first part is highly identical to a regular VMC block.

  -  The second part is to specify the correlated sampling options and
     define the cost function.

  -  The last part is used to specify the options of different optimizers,
     which can be very distinct from one to another.

VMC run for the optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The VMC calculation for the wavefunction optimization has a strict requirement
that ``samples`` or ``samplesperthread`` must be specified because of the optimizer needs for the stored ``samples``.
The input parameters of this part are identical to the VMC method.

Recommendations:

-  Run the inclusive VMC calculation correctly and efficiently because
   this takes a significant amount of time during optimization. For
   example, make sure the derived ``steps`` per block is 1 and use larger ``substeps`` to
   control the correlation between ``samples``.

-  A reasonable starting wavefunction is necessary. A lot of
   optimization fails because of a bad wavefunction starting point. The
   sign of a bad initial wavefunction includes but is not limited to a
   very long equilibration time, low acceptance ratio, and huge
   variance. The first thing to do after a failed optimization is to
   check the information provided by the VMC calculation via
   ``*.scalar.dat files``.

Correlated sampling and cost function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After generating the samples with VMC, the derivatives of the wavefunction with respect to the parameters are computed for proposing a new set of parameters by optimizers.
And later, a correlated sampling calculation is performed to quickly evaluate values of the cost function on the old set of parameters and the new set for further decisions.
The input parameters are listed in the following table.

``linear`` method:

  parameters:

  +--------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | **Name**                 | **Datatype** | **Values**  | **Default** | **Description**                                  |
  +==========================+==============+=============+=============+==================================================+
  | ``nonlocalpp``           | text         | yes, no     | no          | include non-local PP energy in the cost function |
  +--------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``use_nonlocalpp_deriv`` | text         | yes, no     | yes         | Add non-local PP energy derivative contribution  |
  +--------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``minwalkers``           | real         | 0--1        | 0.3         | Lower bound of the effective weight              |
  +--------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``maxWeight``            | real         | :math:`> 1` | 1e6         | Maximum weight allowed in reweighting            |
  +--------------------------+--------------+-------------+-------------+--------------------------------------------------+

Additional information:

- ``maxWeight`` The default should be good.

- ``nonlocalpp`` The ``nonlocalpp`` contribution to the local energy depends on the
  wavefunction. When a new set of parameters is proposed, this
  contribution needs to be updated if the cost function consists of local
  energy. Fortunately, nonlocal contribution is chosen small when making a
  PP for small locality error. We can ignore its change and avoid the
  expensive computational cost. An implementation issue with legacy GPU code is
  that a large amount of memory is consumed with this option.

- ``minwalkers`` This is a ``critical`` parameter. When the ratio of effective samples to actual number of samples in a reweighting step goes lower than ``minwalkers``,
  the proposed set of parameters is invalid.

The cost function consists of three components: energy, unreweighted variance, and reweighted variance.

::

     <cost name="energy">                   0.95 </cost>
     <cost name="unreweightedvariance">     0.00 </cost>
     <cost name="reweightedvariance">       0.05 </cost>

Optimizers
~~~~~~~~~~

QMCPACK implements a number of different optimizers each with different
priorities for accuracy, convergence, memory usage, and stability. The
optimizers can be switched among “OneShiftOnly” (default), “adaptive,”
“descent,” “hybrid,” and “quartic” (old) using the following line in the
optimization block:

::

<parameter name="MinMethod"> THE METHOD YOU LIKE </parameter>

OneShiftOnly Optimizer
~~~~~~~~~~~~~~~~~~~~~~

The OneShiftOnly optimizer targets a fast optimization by moving parameters more aggressively. It works with OpenMP and GPU and can be considered for large systems.
This method relies on the effective weight of correlated sampling rather than the cost function value to justify a new set of parameters.
If the effective weight is larger than ``minwalkers``, the new set is taken whether or not the cost function value decreases.
If a proposed set is rejected, the standard output prints the measured ratio of effective samples to the total number of samples
and adjustment on ``minwalkers`` can be made if needed.

``linear`` method:

  parameters:

  +--------------+--------------+-------------+-------------+---------------------------------------------------+
  | **Name**     | **Datatype** | **Values**  | **Default** | **Description**                                   |
  +==============+==============+=============+=============+===================================================+
  | ``shift_i``  | real         | :math:`> 0` | 0.01        | Direct stabilizer added to the Hamiltonian matrix |
  +--------------+--------------+-------------+-------------+---------------------------------------------------+
  | ``shift_s``  | real         | :math:`> 0` | 1.00        | Initial stabilizer based on the overlap matrix    |
  +--------------+--------------+-------------+-------------+---------------------------------------------------+

Additional information:

-  ``shift_i`` This is the direct term added to the diagonal of the Hamiltonian
   matrix. It provides more stable but slower optimization with a large
   value.

-  ``shift_s`` This is the initial value of the stabilizer based on the overlap
   matrix added to the Hamiltonian matrix. It provides more stable but
   slower optimization with a large value. The used value is
   auto-adjusted by the optimizer.

Recommendations:

- Default ``shift_i``, ``shift_s`` should be fine.

- For hard cases, increasing ``shift_i`` (by a factor of 5 or 10) can significantly stabilize the optimization by reducing the pace towards the optimal parameter set.

- If the VMC energy of the last optimization iterations grows significantly, increase ``minwalkers`` closer to 1 and make the optimization stable.

- If the first iterations of optimization are rejected on a reasonable initial wavefunction,
  lower the ``minwalkers`` value based on the measured value printed in the standard output to accept the move.

We recommended using this optimizer in two sections with a very small ``minwalkers`` in the first and a large value in the second, such as the following.
In the very beginning, parameters are far away from optimal values and large changes are proposed by the optimizer.
Having a small ``minwalkers`` makes it much easier to accept these changes.
When the energy gradually converges, we can have a large ``minwalkers`` to avoid risky parameter sets.

::

  <loop max="6">
   <qmc method="linear" move="pbyp" gpu="yes">
     <!-- Specify the VMC options -->
     <parameter name="walkers">                1 </parameter>
     <parameter name="samples">            10000 </parameter>
     <parameter name="stepsbetweensamples">    1 </parameter>
     <parameter name="substeps">               5 </parameter>
     <parameter name="warmupSteps">            5 </parameter>
     <parameter name="blocks">                25 </parameter>
     <parameter name="timestep">             1.0 </parameter>
     <parameter name="usedrift">              no </parameter>
     <estimator name="LocalEnergy" hdf5="no"/>
     <!-- Specify the optimizer options -->
     <parameter name="MinMethod">    OneShiftOnly </parameter>
     <parameter name="minwalkers">           1e-4 </parameter>
   </qmc>
  </loop>
  <loop max="12">
   <qmc method="linear" move="pbyp" gpu="yes">
     <!-- Specify the VMC options -->
     <parameter name="walkers">                1 </parameter>
     <parameter name="samples">            20000 </parameter>
     <parameter name="stepsbetweensamples">    1 </parameter>
     <parameter name="substeps">               5 </parameter>
     <parameter name="warmupSteps">            2 </parameter>
     <parameter name="blocks">                50 </parameter>
     <parameter name="timestep">             1.0 </parameter>
     <parameter name="usedrift">              no </parameter>
     <estimator name="LocalEnergy" hdf5="no"/>
     <!-- Specify the optimizer options -->
     <parameter name="MinMethod">    OneShiftOnly </parameter>
     <parameter name="minwalkers">            0.5 </parameter>
   </qmc>
  </loop>

For each optimization step, you will see

::

  The new set of parameters is valid. Updating the trial wave function!

or

::

  The new set of parameters is not valid. Revert to the old set!

Occasional rejection is fine. Frequent rejection indicates potential
problems, and users should inspect the VMC calculation or change
optimization strategy. To track the progress of optimization, use the
command ``qmca -q ev *.scalar.dat`` to look at the VMC energy and
variance for each optimization step.

Adaptive Optimizer
~~~~~~~~~~~~~~~~~~

The default setting of the adaptive optimizer is to construct the linear
method Hamiltonian and overlap matrices explicitly and add different
shifts to the Hamiltonian matrix as “stabilizers.” The generalized
eigenvalue problem is solved for each shift to obtain updates to the
wavefunction parameters. Then a correlated sampling is performed for
each shift’s updated wavefunction and the initial trial wavefunction
using the middle shift’s updated wavefunction as the guiding function.
The cost function for these wavefunctions is compared, and the update
corresponding to the best cost function is selected. In the next
iteration, the median magnitude of the stabilizers is set to the
magnitude that generated the best update in the current iteration, thus
adapting the magnitude of the stabilizers automatically.

When the trial wavefunction contains more than 10,000 parameters,
constructing and storing the linear method matrices could become a
memory bottleneck. To avoid explicit construction of these matrices, the
adaptive optimizer implements the block linear method (BLM) approach.
:cite:`Zhao:2017:blocked_lm` The BLM tries to find an
approximate solution :math:`\vec{c}_{opt}` to the standard LM
generalized eigenvalue problem by dividing the variable space into a
number of blocks and making intelligent estimates for which directions
within those blocks will be most important for constructing
:math:`\vec{c}_{opt}`, which is then obtained by solving a smaller, more
memory-efficient eigenproblem in the basis of these supposedly important
block-wise directions.

``linear`` method:

  parameters:

  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | **Name**                  | **Datatype** | **Values**              | **Default** | **Description**                                                                                 |
  +===========================+==============+=========================+=============+=================================================================================================+
  | ``max_relative_change``   | real         | :math:`> 0`             | 10.0        | Allowed change in cost function                                                                 |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``max_param_change``      | real         | :math:`> 0`             | 0.3         | Allowed change in wavefunction parameter                                                        |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``shift_i``               | real         | :math:`> 0`             | 0.01        | Initial diagonal stabilizer added to the Hamiltonian matrix                                     |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``shift_s``               | real         | :math:`> 0`             | 1.00        | Initial overlap-based stabilizer added to the Hamiltonian matrix                                |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``target_shift_i``        | real         | any                     | -1.0        | Diagonal stabilizer value aimed for during adaptive method (disabled if :math:`\leq 0`)         |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``cost_increase_tol``     | real         | :math:`\geq 0`          | 0.0         |  Tolerance for cost function increases                                                          |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``chase_lowest``          | text         | yes, no                 | yes         | Chase the lowest eigenvector in iterative solver                                                |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``chase_closest``         | text         | yes, no                 | no          | Chase the eigenvector closest to initial guess                                                  |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``block_lm``              | text         | yes, no                 | no          | Use BLM                                                                                         |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``blocks``                | integer      | :math:`> 0`             |             | Number of blocks in BLM                                                                         |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``nolds``                 | integer      | :math:`> 0`             |             | Number of old update vectors used in BLM                                                        |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+
  | ``nkept``                 | integer      | :math:`> 0`             |             | Number of eigenvectors to keep per block in BLM                                                 |
  +---------------------------+--------------+-------------------------+-------------+-------------------------------------------------------------------------------------------------+

Additional information:

-  ``shift_i`` This is the initial coefficient used to scale the diagonal
   stabilizer. More stable but slower optimization is expected with a
   large value. The adaptive method will automatically adjust this value
   after each linear method iteration.

-  ``shift_s`` This is the initial coefficient used to scale the overlap-based
   stabilizer. More stable but slower optimization is expected with a
   large value. The adaptive method will automatically adjust this value
   after each linear method iteration.

-  ``target_shift_i`` If set greater than zero, the adaptive method will choose the
   update whose shift_i value is closest to this target value so long as
   the associated cost is within cost_increase_tol of the lowest cost.
   Disable this behavior by setting target_shift_i to a negative number.

-  ``cost_increase_tol`` Tolerance for cost function increases when selecting the best
   shift.

-  ``nblocks`` This is the number of blocks used in BLM. The amount of memory
   required to store LM matrices decreases as the number of blocks
   increases. But the error introduced by BLM would increase as the
   number of blocks increases.

-  ``nolds`` In BLM, the interblock correlation is accounted for by including a
   small number of wavefunction update vectors outside the block. Larger
   would include more interblock correlation and more accurate results
   but also higher memory requirements.

-  ``nkept`` This is the number of update directions retained from each block in
   the BLM. If all directions are retained in each block, then the BLM
   becomes equivalent to the standard LM. Retaining five or fewer
   directions per block is often sufficient.

Recommendations:

-  Default ``shift_i``, ``shift_s`` should be fine.

-  When there are fewer than about 5,000 variables being optimized, the
   traditional LM is preferred because it has a lower overhead than the
   BLM when the number of variables is small.

-  Initial experience with the BLM suggests that a few hundred blocks
   and a handful of and often provide a good balance between memory use
   and accuracy. In general, using fewer blocks should be more accurate
   but would require more memory.

::

  <loop max="15">
   <qmc method="linear" move="pbyp">
     <!-- Specify the VMC options -->
     <parameter name="walkers">                1 </parameter>
     <parameter name="samples">            20000 </parameter>
     <parameter name="stepsbetweensamples">    1 </parameter>
     <parameter name="substeps">               5 </parameter>
     <parameter name="warmupSteps">            5 </parameter>
     <parameter name="blocks">                50 </parameter>
     <parameter name="timestep">             1.0 </parameter>
     <parameter name="usedrift">              no </parameter>
     <estimator name="LocalEnergy" hdf5="no"/>
     <!-- Specify the correlated sampling options and define the cost function -->
          <cost name="energy">               1.00 </cost>
          <cost name="unreweightedvariance"> 0.00 </cost>
          <cost name="reweightedvariance">   0.00 </cost>
     <!-- Specify the optimizer options -->
     <parameter name="MinMethod">adaptive</parameter>
     <parameter name="max_relative_cost_change">10.0</parameter>
     <parameter name="shift_i"> 1.00 </parameter>
     <parameter name="shift_s"> 1.00 </parameter>
     <parameter name="max_param_change"> 0.3 </parameter>
     <parameter name="chase_lowest"> yes </parameter>
     <parameter name="chase_closest"> yes </parameter>
     <parameter name="block_lm"> no </parameter>
     <!-- Specify the BLM specific options if needed
       <parameter name="nblocks"> 100 </parameter>
       <parameter name="nolds"> 5 </parameter>
       <parameter name="nkept"> 3 </parameter>
     -->
   </qmc>
  </loop>

The adaptive optimizer is also able to optimize individual excited states directly. :cite:`Zhao:2016:dir_tar`
In this case, it tries to minimize the following function:

.. math:: \Omega[\Psi]=\frac{\left<\Psi|\omega-H|\Psi\right>}{\left<\Psi|{\left(\omega-H\right)}^2|\Psi\right>}\:.

The global minimum of this function corresponds to the state whose
energy lies immediately above the shift parameter :math:`\omega` in the
energy spectrum. For example, if :math:`\omega` were placed in between
the ground state energy and the first excited state energy and the
wavefunction ansatz was capable of a good description for the first
excited state, then the wavefunction would be optimized for the first
excited state. Note that if the ansatz is not capable of a good
description of the excited state in question, the optimization could
converge to a different state, as is known to occur in some
circumstances for traditional ground state optimizations. Note also that
the ground state can be targeted by this method by choosing
:math:`\omega` to be below the ground state energy, although we should
stress that this is not the same thing as a traditional ground state
optimization and will in general give a slightly different wavefunction.
Excited state targeting requires two additional parameters, as shown in
the following table.

Excited state targeting:

  parameters:

  +-------------------+--------------+--------------+-------------+---------------------------------------------------------+
  | **Name**          | **Datatype** | **Values**   | **Default** | **Description**                                         |
  +===================+==============+==============+=============+=========================================================+
  | ``targetExcited`` | text         | yes, no      | no          | Whether to use the excited state targeting optimization |
  +-------------------+--------------+--------------+-------------+---------------------------------------------------------+
  | ``omega``         | real         | real numbers | none        | Energy shift used to target different excited states    |
  +-------------------+--------------+--------------+-------------+---------------------------------------------------------+

Excited state recommendations:

-  Because of the finite variance in any approximate wavefunction, we
   recommended setting :math:`\omega=\omega_0-\sigma`, where
   :math:`\omega_0` is placed just below the energy of the targeted
   state and :math:`\sigma^2` is the energy variance.

-  To obtain an unbiased excitation energy, the ground state should be
   optimized with the excited state variational principle as well by
   setting ``omega`` below the ground state energy. Note that using the ground
   state variational principle for the ground state and the excited
   state variational principle for the excited state creates a bias in
   favor of the ground state.

Descent Optimizer
~~~~~~~~~~~~~~~~~

Gradient descent algorithms are an alternative set of optimization methods to the OneShiftOnly and adaptive optimizers based on the linear method.
These methods use only first derivatives to optimize trial wave functions and convergence can be accelerated by retaining a memory of previous derivative values.
Multiple flavors of accelerated descent methods are available. They differ in details such as the schemes for adaptive adjustment of step sizes. :cite:`Otis2019`
Descent algorithms avoid the construction of matrices that occurs in the linear method and consequently can be applied to larger sets of
optimizable parameters.
Parameters for descent are shown in the table below.

``descent`` method:

  parameters:

  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | **Name**            | **Datatype** | **Values**                     | **Default** | **Description**                                                 |
  +=====================+==============+================================+=============+=================================================================+
  | ``flavor``          | text         | RMSprop, Random, ADAM, AMSGrad | RMSprop     | Particular type of descent method                               |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``Ramp_eta``        | text         | yes, no                        | no          | Whether to gradually ramp up step sizes                         |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``Ramp_num``        | integer      | :math:`> 0`                    | 30          | Number of steps over which to ramp up step size                 |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``TJF_2Body_eta``   | real         | :math:`> 0`                    | 0.01        | Step size for two body Jastrow parameters                       |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``TJF_1Body_eta``   | real         | :math:`> 0`                    | 0.01        | Step size for one body Jastrow parameters                       |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``F_eta``           | real         | :math:`> 0`                    | 0.001       | Step size for number counting Jastrow F matrix parameters       |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``Gauss_eta``       | real         | :math:`> 0`                    | 0.001       | Step size for number counting Jastrow gaussian basis parameters |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``CI_eta``          | real         | :math:`> 0`                    | 0.01        | Step size for CI parameters                                     |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``Orb_eta``         | real         | :math:`> 0`                    | 0.001       | Step size for orbital parameters                                |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``collection_step`` | real         | :math:`> 0`                    | 0.01        | Step number to start collecting samples for final averages      |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``compute_step``    | real         | :math:`> 0`                    | 0.001       | Step number to start computing averaged from stored history     |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+
  | ``print_derivs``    | real         | yes, no                        | no          | Whether to print parameter derivatives                          |
  +---------------------+--------------+--------------------------------+-------------+-----------------------------------------------------------------+


These descent algorithms have been extended to the optimization of the same excited state functional as the adaptive LM. :cite:`Otis2020`
This also allows the hybrid optimizer discussed below to be applied to excited states.
The relevant parameters are the same as for targeting excited states with the adaptive optimizer above.

Additional information and recommendations:

-  It is generally advantageous to set different step sizes for
   different types of parameters. More nonlinear parameters such as
   those for number counting Jastrow factors or orbitals typically
   require smaller steps sizes than those for CI coefficients or
   traditional Jastrow parameters. There are defaults for several
   parameter types and a default of .001 has been chosen for all other
   parameters.

-  The ability to gradually ramp up step sizes to their input values is
   useful for avoiding spikes in the average local energy during early
   iterations of descent optimization. This initial rise in the energy
   occurs as a memory of past gradients is being built up and it may be
   possible for the energy to recover without ramping if there are
   enough iterations in the optimization.

-  The step sizes chosen can have a substantial influence on the quality
   of the optimization and the final variational energy achieved. Larger
   step sizes may be helpful if there is reason to think the descent
   optimization is not reaching the minimum energy. There are also
   additional hyperparameters in the descent algorithms with default
   values. :cite:`Otis2019` They seem to have limited
   influence on the effectiveness of the optimization compared to step
   sizes, but users can adjust them within the source code of the
   descent engine if they wish.

-  The sampling effort for individual descent steps can be small
   compared that for linear method iterations as shown in the example
   input below. Something in the range of 10,000 to 30,000 seems
   sufficient for molecules with tens of electrons. However, descent
   optimizations may require anywhere from a few hundred to a few
   thousand iterations.
 
 -  For reporting quantities such as a final energy and associated uncertainty,
    an average over many descent steps can be taken. The parameters for 
    ``collection_step`` and ``compute_step`` help automate this task.
    After the descent iteration specified by ``collection_step``, a 
    history of local energy values will be kept for determining a final 
    error and average, which will be computed and given in the output 
    once the iteration specified by ``compute_step`` is reached. For 
    reasonable results, this procedure should use descent steps near 
    the end of the optimization when the wave function parameters are essentially 
    no longer changing.

-  In cases where a descent optimization struggles to reach the minimum
   and a linear method optimization is not possible or unsatisfactory,
   it may be useful to try the hybrid optimization approach described in
   the next subsection.

::


  <loop max="2000">
     <qmc method="linear" move="pbyp" checkpoint="-1" gpu="no">

     <!-- VMC inputs -->
      <parameter name="blocks">2000</parameter>
      <parameter name="steps">1</parameter>
      <parameter name="samples">20000</parameter>
      <parameter name="warmupsteps">100</parameter>
      <parameter name="timestep">0.05</parameter>

      <parameter name="MinMethod">descent</parameter>
      <estimator name="LocalEnergy" hdf5="no"/>
      <parameter name="usebuffer">yes</parameter>

      <estimator name="LocalEnergy" hdf5="no"/>

      <!-- Descent Inputs -->
        <parameter name="flavor">RMSprop</parameter>

        <parameter name="Ramp_eta">no</parameter>
        <parameter name="Ramp_num">30</parameter>

       <parameter name="TJF_2Body_eta">.02</parameter>
        <parameter name="TJF_1Body_eta">.02</parameter>
       <parameter name="F_eta">.001</parameter>
       <parameter name="Gauss_eta">.001</parameter>
       <parameter name="CI_eta">.1</parameter>
       <parameter name="Orb_eta">.0001</parameter>

       <parameter name="collection_step">500</parameter>
       <parameter name="compute_step">998</parameter>
       
      <parameter name="targetExcited"> yes </parameter>
      <parameter name="targetExcited"> -11.4 </parameter>

       <parameter name="print_derivs">no</parameter>


     </qmc>
  </loop>

Hybrid Optimizer
~~~~~~~~~~~~~~~~

Another optimization option is to use a hybrid combination of accelerated descent and blocked linear method.
It provides a means to retain the advantages of both individual methods while scaling to large numbers of parameters beyond the traditional 10,000 parameter limit of the linear method. :cite:`Otis2019`
In a hybrid optimization, alternating sections of descent and BLM optimization are used.
Gradient descent is used to identify the previous important directions in parameter space used by the BLM, the number of which is set by the ``nold`` input for the BLM.
Over the course of a section of descent, vectors of parameter differences are stored and then passed to the linear method engine after the optimization changes to the BLM.
One motivation for including sections of descent is to counteract noise in linear method updates due to uncertainties in its step direction and allow for a smoother movement to the minimum.
There are two additional parameters used in the hybrid optimization and it requires a slightly different format of input to specify the constituent methods as shown below in the example.

``descent`` method:

  parameters:

  +---------------------+--------------+-------------+-------------+--------------------------------------+
  | **Name**            | **Datatype** | **Values**  | **Default** | **Description**                      |
  +=====================+==============+=============+=============+======================================+
  | ``num_updates``     | integer      | :math:`> 0` |             | Number of steps for a method         |
  +---------------------+--------------+-------------+-------------+--------------------------------------+
  | ``Stored_Vectors``  | integer      | :math:`> 0` | 5           | Number of vectors to transfer to BLM |
  +---------------------+--------------+-------------+-------------+--------------------------------------+

::


  <loop max="203">
  <qmc method="linear" move="pbyp" checkpoint="-1" gpu="no">
   <parameter name="Minmethod"> hybrid </parameter>

   <optimizer num_updates="100">

  <parameter name="blocks">1000</parameter>
       <parameter name="steps">1</parameter>
       <parameter name="samples">20000</parameter>
       <parameter name="warmupsteps">1000</parameter>
       <parameter name="timestep">0.05</parameter>

       <estimator name="LocalEnergy" hdf5="no"/>

       <parameter name="Minmethod"> descent </parameter>
       <parameter name="Stored_Vectors">5</parameter>
       <parameter name="flavor">RMSprop</parameter>
       <parameter name="TJF_2Body_eta">.01</parameter>
       <parameter name="TJF_1Body_eta">.01</parameter>
       <parameter name="CI_eta">.1</parameter>

       <parameter name="Ramp_eta">no</parameter>
       <parameter name="Ramp_num">10</parameter>
   </optimizer>

   <optimizer num_updates="3">

       <parameter name="blocks">2000</parameter>
       <parameter name="steps">1</parameter>
       <parameter name="samples">1000000</parameter>
       <parameter name="warmupsteps">1000</parameter>
       <parameter name="timestep">0.05</parameter>

       <estimator name="LocalEnergy" hdf5="no"/>

       <parameter name="Minmethod"> adaptive </parameter>
       <parameter name="max_relative_cost_change">10.0</parameter>
       <parameter name="max_param_change">3</parameter>
       <parameter name="shift_i">0.01</parameter>
       <parameter name="shift_s">1.00</parameter>

       <parameter name="block_lm">yes</parameter>
       <parameter name="nblocks">2</parameter>
       <parameter name="nolds">5</parameter>
       <parameter name="nkept">5</parameter>

   </optimizer>
  </qmc>
  </loop>

Additional information and recommendations:

-  In the example above, the input for ``loop`` gives the total number
   of steps for the full optimization while the inputs for
   ``num_updates`` specify the number of steps in the constituent
   methods. For this case, the optimization would begin with 100 steps
   of descent using the parameters in the first ``optimizer`` block and
   then switch to the BLM for 3 steps before switching back to descent
   for the final 100 iterations of the total of 203.

-  The design of the hybrid method allows for more than two
   ``optimizer`` blocks to be used and the optimization will cycle
   through the individual methods. However, the effectiveness of this in
   terms of the quality of optimization results is unexplored.

-  It can be useful to follow a hybrid optimization with a section of
   pure descent optimization and take an average energy over the last
   few hundred iterations as the final variational energy. This approach
   can achieve a lower statistical uncertainty on the energy for less
   overall sampling effort compared to what a pure linear method
   optimization would require. The ``collection_step`` and ``compute_step``
   parameters discussed earlier for descent are useful for setting up
   the descent engine to do this averaging on its own.

Quartic Optimizer
~~~~~~~~~~~~~~~~~

*This is an older optimizer method retained for compatibility. We
recommend starting with the newest OneShiftOnly or adaptive optimizers.*
The quartic optimizer fits a quartic polynomial to 7 values of the cost
function obtained using reweighting along the chosen direction and
determines the optimal move. This optimizer is very robust but is a bit
conservative when accepting new steps, especially when large parameters
changes are proposed.

``linear`` method:

  parameters:

  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+
  | **Name**              | **Datatype** | **Values**  | **Default** | **Description**                                  |
  +=======================+==============+=============+=============+==================================================+
  | ``bigchange``         | real         | :math:`> 0` | 50.0        | Largest parameter change allowed                 |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``alloweddifference`` | real         | :math:`> 0` | 1e-4        | Allowed increase in energy                       |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``exp0``              | real         | any value   | -16.0       | Initial value for stabilizer                     |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``stabilizerscale``   | real         | :math:`> 0` | 2.0         | Increase in value of ``exp0`` between iterations |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``nstabilizers``      | integer      | :math:`> 0` | 3           | Number of stabilizers to try                     |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``max_its``           | integer      | :math:`> 0` | 1           | Number of inner loops with same samples          |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------------+

Additional information:

-  ``exp0`` This is the initial value for stabilizer (shift to diagonal of H).
   The actual value of stabilizer is :math:`10^{\textrm{exp0}}`.

Recommendations:

-  For hard cases (e.g., simultaneous optimization of long MSD and
   3-Body J), set ``exp0`` to 0 and do a single inner iteration (max its=1) per
   sample of configurations.

::

  <!-- Specify the optimizer options -->
  <parameter name="MinMethod">quartic</parameter>
  <parameter name="exp0">-6</parameter>
  <parameter name="alloweddifference"> 1.0e-4 </parameter>
  <parameter name="nstabilizers"> 1 </parameter>
  <parameter name="bigchange">15.0</parameter>

General Recommendations
~~~~~~~~~~~~~~~~~~~~~~~

-  All electron wavefunctions are typically more difficult to optimize
   than pseudopotential wavefunctions because of the importance of the
   wavefunction near the nucleus.

-  Two-body Jastrow contributes the largest portion of correlation
   energy from bare Slater determinants. Consequently, the recommended
   order for optimizing wavefunction components is two-body, one-body,
   three-body Jastrow factors and MSD coefficients.

-  For two-body spline Jastrows, always start from a reasonable one. The
   lack of physically motivated constraints in the functional form at
   large distances can cause slow convergence if starting from zero.

-  One-body spline Jastrow from old calculations can be a good starting
   point.

-  Three-body polynomial Jastrow can start from zero. It is beneficial
   to first optimize one-body and two-body Jastrow factors without
   adding three-body terms in the calculation and then add the
   three-body Jastrow and optimize all the three components together.

Optimization of CI coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When storing a CI wavefunction in HDF5 format, the CI coefficients and
the :math:`\alpha` and :math:`\beta` components of each CI are not in
the XML input file. When optimizing the CI coefficients, they will be
stored in HDF5 format. The optimization header block will have to
specify that the new CI coefficients will be saved to HDF5 format. If
the tag is not added coefficients will not be saved.

::

  <qmc method="linear" move="pbyp" gpu="no" hdf5="yes">

  The rest of the optimization block remains the same.

When running the optimization, the new coefficients will be stored in a ``*.sXXX.opt.h5`` file,  where XXX coressponds to the series number. The H5 file contains only the optimized coefficients. The corresponding ``*.sXXX.opt.xml`` will be updated for each optimization block as follows:

::

  <detlist size="1487" type="DETS" nca="0" ncb="0" nea="2" neb="2" nstates="85" cutoff="1e-2" href="../LiH.orbs.h5" opt_coeffs="LiH.s001.opt.h5"/>

The opt_coeffs tag will then reference where the new CI coefficients are
stored.

When restarting the run with the new optimized coeffs, you need to
specify the previous hdf5 containing the basis set, orbitals, and MSD,
as well as the new optimized coefficients. The code will read the
previous data but will rewrite the coefficients that were optimized with
the values found in the \*.sXXX.opt.h5 file. Be careful to keep the pair
of optimized CI coefficients and Jastrow coefficients together to avoid
inconsistencies.

Parameter gradients
~~~~~~~~~~~~~~~~~~~
The gradients of the energy with respect to the variational parameters can be checked and optionally written to a file.
The check compares the analytic derivatives with a finite difference approximation.
These are activated by giving a ``gradient_test`` method in and ``optimize`` block, as follows:

::

     <qmc method="linear" move="pbyp">
      <optimize method="gradient_test">
      </optimize>
      ... rest of optimizer input ...

The check will print a table to the standard output with the parameter name, value, analytic gradient, finite difference gradient, and the percent difference between them.

Writing the analytic parameter gradients to a file is enabled by using the ``output_param_file`` parameter.
The file name is ``<project id>.param.s000.scalar.dat``.
It contains one line per loop iteration, to allow using existing tools to compute averages and error bars on the values.

  +-----------------------+--------------+-------------+-------------+--------------------------------------------+
  | **Name**              | **Datatype** | **Values**  | **Default** | **Description**                            |
  +=======================+==============+=============+=============+============================================+
  | ``output_param_file`` | text         | yes, no     | no          |  Output parameter gradients to a file      |
  +-----------------------+--------------+-------------+-------------+--------------------------------------------+

The input would look like the following:

::

    <qmc method="linear" move="pbyp" checkpoint="-1" gpu="no">
      <optimize method="gradient_test">
        <parameter name="output_param_file">yes</parameter>
      </optimize>
      ... rest of optimizer input ...



Output of intermediate values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following parameters to the linear optimizers to output intermediate values such as the overlap and Hamiltonian matrices.

  +-------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | **Name**                | **Datatype** | **Values**  | **Default** | **Description**                                  |
  +=========================+==============+=============+=============+==================================================+
  | ``output_matrices_csv`` | text         | yes, no     | no          |  Output linear method matrices to CSV files      |
  +-------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``output_matrices_hdf`` | text         | yes, no     | no          |  Output linear method matrices to HDF file       |
  +-------------------------+--------------+-------------+-------------+--------------------------------------------------+
  | ``freeze_parameters``   | text         | yes, no     | no          |  Do not update parameters between iterations     |
  +-------------------------+--------------+-------------+-------------+--------------------------------------------------+

  The ``output_matrices_csv`` parameter will write to <base name>.ham.s000.scalar.dat and <base name>.ovl.scalar.dat.  One line per iteration of the optimizer loop.  Combined with ``freeze_parameters``, this allows computing error bars on the matrices for use in regression testing.

  The ``output_matrices_hdf`` parameter will output in HDF format the matrices used in the linear method along with the shifts and the eigenvalue and eigenvector produced by QMCPACK.  The file is named "<base name>.<series number>.linear_matrices.h5".  It only works with the batched optimizer (batched version of ``linear``)


.. _dmc:

Diffusion Monte Carlo
---------------------

``dmc`` driver
~~~~~~~~~~~~~~

Main input parameters are given in :numref:`table9`, additional in :numref:`table10`.

parameters:

.. _table9:
.. table::

  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | **Name**                       | **Datatype** | **Values**              | **Default** | **Description**                               |
  +================================+==============+=========================+=============+===============================================+
  | ``targetwalkers``              | integer      | :math:`> 0`             | dep.        | Overall total number of walkers               |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``blocks``                     | integer      | :math:`\geq 0`          | 1           | Number of blocks                              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``steps``                      | integer      | :math:`\geq 0`          | 1           | Number of steps per block                     |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``warmupsteps``                | integer      | :math:`\geq 0`          | 0           | Number of steps for warming up                |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``timestep``                   | real         | :math:`> 0`             | 0.1         | Time step for each electron move              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``nonlocalmoves``              | string       | yes, no, v0, v1, v3     | no          | Run with T-moves                              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``branching_cutoff_scheme``    |              |                         |             |                                               |
  |                                |              |                         |             |                                               |
  |                                | string       | classic/DRV/ZSGMA/YL    | classic     | Branch cutoff scheme                          |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``maxcpusecs``                 | real         | :math:`\geq 0`          | 3.6e5       | Deprecated. Superseded by ``max_seconds``     |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``max_seconds``                | real         | :math:`\geq 0`          | 3.6e5       | Maximum allowed walltime in seconds           |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``blocks_between_recompute``   | integer      | :math:`\geq 0`          | dep.        | Wavefunction recompute frequency              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``spinMass``                   | real         | :math:`> 0`             | 1.0         | Effective mass for spin sampling              |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+
  | ``debug_checks``               | text         | see additional info     | dep.        | Turn on/off additional recompute and checks   |
  +--------------------------------+--------------+-------------------------+-------------+-----------------------------------------------+

.. centered:: Table 9 Main DMC input parameters.

.. _table10:
.. table::

  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | **Name**                    | **Datatype** | **Values**              | **Default** | **Description**                         |
  +=============================+==============+=========================+=============+=========================================+
  | ``energyUpdateInterval``    | integer      | :math:`\geq 0`          | 0           | Trial energy update interval            |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``refEnergy``               | real         | all values              | dep.        | Reference energy in atomic units        |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``feedback``                | double       | :math:`\geq 0`          | 1.0         | Population feedback on the trial energy |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``sigmaBound``              | 10           | :math:`\geq 0`          | 10          | Parameter to cutoff large weights       |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``killnode``                | string       | yes/other               | no          | Kill or reject walkers that cross nodes |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``warmupByReconfiguration`` | option       | yes,no                  | 0           | Warm up with a fixed population         |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``reconfiguration``         | string       | yes/pure/other          | no          | Fixed population technique              |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``branchInterval``          | integer      | :math:`\geq 0`          | 1           | Branching interval                      |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``substeps``                | integer      | :math:`\geq 0`          | 1           | Branching interval                      |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``MaxAge``                  | double       | :math:`\geq 0`          | 10          | Kill persistent walkers                 |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``MaxCopy``                 | double       | :math:`\geq 0`          | 2           | Limit population growth                 |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``maxDisplSq``              | real         | all values              | -1          | Maximum particle move                   |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``scaleweight``             | string       | yes/other               | yes         | Scale weights (CUDA only)               |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``checkproperties``         | integer      | :math:`\geq 0`          | 100         | Number of steps between walker updates  |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``fastgrad``                | text         | yes/other               | yes         | Fast gradients                          |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``storeconfigs``            | integer      | all values              | 0           | Store configurations                    |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``use_nonblocking``         | string       | yes/no                  | yes         | Using nonblocking send/recv             |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+
  | ``debug_disable_branching`` | string       | yes/no                  | no          | Disable branching for debugging         |
  |                             |              |                         |             | without correctness guarantee           |
  +-----------------------------+--------------+-------------------------+-------------+-----------------------------------------+

.. centered:: Table 10 Additional DMC input parameters.

Additional information:

-  ``targetwalkers``: A DMC run can be considered a restart run or a new
   run. A restart run is considered to be any method block beyond the
   first one, such as when a DMC method block follows a VMC block.
   Alternatively, a user reading in configurations from disk would also
   considered a restart run. In the case of a restart run, the DMC
   driver will use the configurations from the previous run, and this
   variable will not be used. For a new run, if the number of walkers is
   less than the number of threads, then the number of walkers will be
   set equal to the number of threads.

-  ``blocks``: This is the number of blocks run during a DMC method
   block. A block consists of a number of DMC steps (steps), after which
   all the statistics accumulated in the block are written to disk.

-  ``steps``: This is the number of DMC steps in a block.

-  ``warmupsteps``: These are the steps at the beginning of a DMC run in
   which the instantaneous average energy is used to update the trial
   energy. During regular steps, E\ :math:`_{ref}` is used.

-  ``timestep``: The ``timestep`` determines the accuracy of the
   imaginary time propagator. Generally, multiple time steps are used to
   extrapolate to the infinite time step limit. A good range of time
   steps in which to perform time step extrapolation will typically have
   a minimum of 99% acceptance probability for each step.

-  ``checkproperties``: When using a particle-by-particle driver, this
   variable specifies how often to reset all the variables kept in the
   buffer.

-  ``maxcpusecs``: Deprecated. Superseded by ``max_seconds``.

-  ``max_seconds``: The default is 100 hours. Once the specified time has
   elapsed, the program will finalize the simulation even if all blocks
   are not completed.

-  ``spinMass`` This is an optional parameter to allow the user to change the rate of spin sampling. If spin sampling is on using ``spinor`` == yes in the electron ParticleSet input, the spin mass determines the rate 
   of spin sampling, resulting in an effective spin timestep :math:`\tau_s = \frac{\tau}{\mu_s}` where 
   :math:`\tau` is the normal spatial timestep and :math:`\mu_s` is the value of the spin mass. The algorithm is described in detail in :cite:`Melton2016-1` and :cite:`Melton2016-2`.

- ``debug_checks`` valid values are 'no', 'all', 'checkGL_after_moves'. If the build type is `debug`, the default value is 'all'. Otherwise, the default value is 'no'.

-  ``energyUpdateInterval``: The default is to update the trial energy
   at every step. Otherwise the trial energy is updated every
   ``energyUpdateInterval`` step.

.. math::

  E_{\text{trial}}=
  \textrm{refEnergy}+\textrm{feedback}\cdot(\ln\texttt{targetWalkers}-\ln N)\:,

where :math:`N` is the current population.

-  ``refEnergy``: The default reference energy is taken from the VMC run
   that precedes the DMC run. This value is updated to the current mean
   whenever branching happens.

-  ``feedback``: This variable is used to determine how strong to react
   to population fluctuations when doing population control. See the
   equation in energyUpdateInterval for more details.

-  ``useBareTau``: The same time step is used whether or not a move is
   rejected. The default is to use an effective time step when a move is
   rejected.

-  ``warmupByReconfiguration``: Warmup DMC is done with a fixed
   population.

-  ``sigmaBound``: This determines the branch cutoff to limit wild
   weights based on the sigma and ``sigmaBound``.

-  ``killnode``: When running fixed-node, if a walker attempts to cross
   a node, the move will normally be rejected. If ``killnode`` = “yes,"
   then walkers are destroyed when they cross a node.

-  ``reconfiguration``: If ``reconfiguration`` is “yes," then run with a
   fixed walker population using the reconfiguration technique.

-  ``branchInterval``: This is the number of steps between branching.
   The total number of DMC steps in a block will be
   ``BranchInterval``\ \*Steps.

-  ``substeps``: This is the same as ``BranchInterval``.

-  ``nonlocalmoves``: Evaluate pseudopotentials using one of the
   nonlocal move algorithms such as T-moves.

   -  no(default): Imposes the locality approximation.

   -  yes/v0: Implements the algorithm in the 2006 Casula
      paper :cite:`Casula2006`.

   -  v1: Implements the v1 algorithm in the 2010 Casula
      paper :cite:`Casula2010`.

   -  v2: Is **not implemented** and is **skipped** to avoid any confusion
      with the v2 algorithm in the 2010 Casula
      paper :cite:`Casula2010`.

   -  v3: (Experimental) Implements an algorithm similar to v1 but is much
      faster. v1 computes the transition probability before each single
      electron T-move selection because of the acceptance of previous
      T-moves. v3 mostly reuses the transition probability computed during
      the evaluation of nonlocal pseudopotentials for the local energy,
      namely before accepting any T-moves, and only recomputes the
      transition probability of the electrons within the same
      pseudopotential region of any electrons touched by T-moves. This is
      an approximation to v1 and results in a slightly different time step
      error, but it significantly reduces the computational cost. v1 and v3
      agree at zero time step. This faster algorithm is the topic of a
      paper in preparation.

      The v1 and v3 algorithms are size-consistent and are important advances over the previous v0 non-size-consistent algorithm. We highly recommend investigating the importance of size-consistency.


-  ``scaleweight``: This is the scaling weight per Umrigar/Nightingale.
   CUDA only.

-  ``MaxAge``: Set the weight of a walker to min(currentweight,0.5)
   after a walker has not moved for ``MaxAge`` steps. Needed if
   persistent walkers appear during the course of a run.

-  ``MaxCopy``: When determining the number of copies of a walker to
   branch, set the number of copies equal to min(Multiplicity,MaxCopy).

-  ``fastgrad``: This calculates gradients with either the fast version
   or the full-ratio version.

-  ``maxDisplSq``: When running a DMC calculation with particle by
   particle, this sets the maximum displacement allowed for a single
   particle move. All distance displacements larger than the max are
   rejected. If initialized to a negative value, it becomes equal to
   Lattice(LR/rc).

-  ``sigmaBound``: This determines the branch cutoff to limit wild
   weights based on the sigma and ``sigmaBound``.

-  ``storeconfigs``: If ``storeconfigs`` is set to a nonzero value, then
   electron configurations during the DMC run will be saved. This option
   is disabled for the OpenMP version of DMC.

-  ``blocks_between_recompute``: See details in :ref:`vmc`.

-  ``branching_cutoff_scheme:`` Modifies how the branching factor is
   computed so as to avoid divergences and stability problems near nodal
   surfaces.

   -  classic (default): The implementation found in QMCPACK v3.0.0 and
      earlier.
      :math:`E_{\rm cut}=\mathrm{min}(\mathrm{max}(\sigma^2 \times \mathrm{sigmaBound},\mathrm{maxSigma}),2.5/\tau)`,
      where :math:`\sigma^2` is the variance and
      :math:`\mathrm{maxSigma}` is set to 50 during warmup
      (equilibration) and 10 thereafter. :math:`\mathrm{sigmaBound}` is
      default to 10.

   -  DRV: Implements the algorithm of DePasquale et al., Eq. 3 in
      :cite:`DePasqualeReliable1988` or Eq. 9 of
      :cite:`Umrigar1993`.
      :math:`E_{\rm cut}=2.0/\sqrt{\tau}`.

   -  ZSGMA: Implements the “ZSGMA” algorithm of
      :cite:`ZenBoosting2016` with :math:`\alpha=0.2`.
      The cutoff energy is modified by a factor including the electron
      count, :math:`E_{\rm cut}=\alpha \sqrt{N/\tau}`, which greatly
      improves size consistency over Eq. 39 of
      :cite:`Umrigar1993`. See Eq. 6 in
      :cite:`ZenBoosting2016` and for an application to
      molecular crystals :cite:`ZenFast2018`.

   -  YL: An unpublished algorithm due to Ye Luo.
      :math:`E_{\rm cut}=\sigma\times\mathrm{min}(\mathrm{sigmaBound},\sqrt{1/\tau})`.
      This option takes into account both size consistency and
      wavefunction quality via the term :math:`\sigma`.
      :math:`\mathrm{sigmaBound}` is default to 10.

.. code-block::
  :caption: The following is an example of a very simple DMC section.
  :name: Listing 44

  <qmc method="dmc" move="pbyp" target="e">
    <parameter name="blocks">100</parameter>
    <parameter name="steps">400</parameter>
    <parameter name="timestep">0.010</parameter>
    <parameter name="warmupsteps">100</parameter>
  </qmc>

The time step should be individually adjusted for each problem.  Please refer to the theory section
on diffusion Monte Carlo.

.. code-block::
  :caption: The following is an example of running a simulation that can be restarted.
  :name: Listing 45

  <qmc method="dmc" move="pbyp"  checkpoint="0">
    <parameter name="timestep">         0.004  </parameter>
    <parameter name="blocks">           100   </parameter>
    <parameter name="steps">            400    </parameter>
  </qmc>

The checkpoint flag instructs QMCPACK to output walker configurations.
This also works in VMC. This will output an h5 file with the name
``projectid.run-number.config.h5``. Check that this file exists before
attempting a restart. To read in this file for a continuation run,
specify the following:

.. code-block::
  :caption: Restart (read walkers from previous run).
  :name: Listing 46

  <mcwalkerset fileroot="BH.s002" version="0 6" collected="yes"/>

BH is the project id, and s002 is the calculation number to read in the walkers from the previous run.

Combining VMC and DMC in a single run (wavefunction optimization can be combined in this way too) is the standard way in which QMCPACK is typically run.   There is no need to run two separate jobs since method sections can be stacked and walkers are transferred between them.

.. code-block::
  :caption: Combined VMC and DMC run.
  :name: Listing 47

  <qmc method="vmc" move="pbyp" target="e">
    <parameter name="blocks">100</parameter>
    <parameter name="steps">4000</parameter>
    <parameter name="warmupsteps">100</parameter>
    <parameter name="samples">1920</parameter>
    <parameter name="walkers">1</parameter>
    <parameter name="timestep">0.5</parameter>
  </qmc>
  <qmc method="dmc" move="pbyp" target="e">
    <parameter name="blocks">100</parameter>
    <parameter name="steps">400</parameter>
    <parameter name="timestep">0.010</parameter>
    <parameter name="warmupsteps">100</parameter>
  </qmc>
  <qmc method="dmc" move="pbyp" target="e">
    <parameter name="warmupsteps">500</parameter>
    <parameter name="blocks">50</parameter>
    <parameter name="steps">100</parameter>
    <parameter name="timestep">0.005</parameter>
  </qmc>

.. _dmc_batch:

Batched ``dmc`` driver (experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  parameters:

  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | **Name**                       | **Datatype** | **Values**              | **Default** | **Description**                                 |
  +================================+==============+=========================+=============+=================================================+
  | ``total_walkers``              | integer      | :math:`> 0`             | 1           | Total number of walkers over all MPI ranks      |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``walkers_per_rank``           | integer      | :math:`> 0`             | 1           | Number of walkers per MPI rank                  |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``crowds``                     | integer      | :math:`> 0`             | dep.        | Number of desynchronized dwalker crowds         |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``blocks``                     | integer      | :math:`\geq 0`          | 1           | Number of blocks                                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``steps``                      | integer      | :math:`\geq 0`          | 1           | Number of steps per block                       |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``warmupsteps``                | integer      | :math:`\geq 0`          | 0           | Number of steps for warming up                  |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``timestep``                   | real         | :math:`> 0`             | 0.1         | Time step for each electron move                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``nonlocalmoves``              | string       | yes, no, v0, v1, v3     | no          | Run with T-moves                                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``branching_cutoff_scheme``    | string       | classic/DRV/ZSGMA/YL    | classic     | Branch cutoff scheme                            |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``blocks_between_recompute``   | integer      | :math:`\geq 0`          | dep.        | Wavefunction recompute frequency                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``feedback``                   | double       | :math:`\geq 0`          | 1.0         | Population feedback on the trial energy         |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``sigmaBound``                 | 10           | :math:`\geq 0`          | 10          | Parameter to cutoff large weights               |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``reconfiguration``            | string       | yes/pure/other          | no          | Fixed population technique                      |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``storeconfigs``               | integer      | all values              | 0           | Store configurations                            |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``use_nonblocking``            | string       | yes/no                  | yes         | Using nonblocking send/recv                     |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``debug_disable_branching``    | string       | yes/no                  | no          | Disable branching for debugging                 |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``crowd_serialize_walkers``    | integer      | yes, no                 | no          | Force use of single walker APIs (for testing)   |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``debug_checks``               | text         | see additional info     | dep.        | Turn on/off additional recompute and checks     |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``spin_mass``                  | real         | :math:`\geq 0`          | 1.0         | Effective mass for spin sampling                |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+
  | ``measure_imbalance``          | text         | yes,no                  | no          | Measure load imbalance at the end of each block |
  +--------------------------------+--------------+-------------------------+-------------+-------------------------------------------------+


- ``crowds`` The number of crowds that the walkers are subdivided into on each MPI rank. If not provided, it is set equal to the number of OpenMP threads.

- ``walkers_per_rank`` The number of walkers per MPI rank. This number does not have to be a multiple of the number of OpenMP
  threads. However, to avoid any idle resources, it is recommended to be at least the number of OpenMP threads for pure CPU runs.
  For GPU runs, a scan of this parameter is necessary to reach reasonable single rank efficiency and also get a balanced time to
  solution. For highest throughput on GPUs, expect to use hundreds of walkers_per_rank, or the largest number that will fit in GPU
  memory. If neither ``total_walkers`` nor ``walkers_per_rank`` is provided, ``walkers_per_rank`` is set equal to ``crowds``.

- ``total_walkers`` Total number of walkers summed over all MPI ranks, or equivalently the total number of walkers in the DMC
  calculation. If not provided, it is computed as ``walkers_per_rank`` times the number of MPI ranks. If both ``total_walkers``
  and ``walkers_per_rank`` are provided, which is not recommended, ``total_walkers`` must be consistently set equal to
  ``walkers_per_rank`` times the number MPI ranks.

- ``debug_checks`` valid values are 'no', 'all', 'checkGL_after_load', 'checkGL_after_moves', 'checkGL_after_tmove'. If the build type is `debug`, the default value is 'all'. Otherwise, the default value is 'no'.

- ``spin_mass`` Optional parameter to allow the user to change the rate of spin sampling. If spin sampling is on using ``spinor`` == yes in the electron ParticleSet input,  the spin mass determines the rate
  of spin sampling, resulting in an effective spin timestep :math:`\tau_s = \frac{\tau}{\mu_s}`. The algorithm is described in detail in :cite:`Melton2016-1` and :cite:`Melton2016-2`.

.. code-block::
  :caption: The following is an example of a minimal DMC section using the batched ``dmc`` driver
  :name: Listing 48b

  <qmc method="dmc" move="pbyp" target="e">
    <parameter name="walkers_per_rank">256</parameter>
    <parameter name="blocks">100</parameter>
    <parameter name="steps">400</parameter>
    <parameter name="timestep">0.010</parameter>
    <parameter name="warmupsteps">100</parameter>
  </qmc>

.. _rmc:

Reptation Monte Carlo
---------------------

Like DMC, RMC is a projector-based method that allows sampling of the
fixed-node wavefunciton. However, by exploiting the path-integral
formulation of Schrödinger’s equation, the RMC algorithm can offer some
advantages over traditional DMC, such as sampling both the mixed and
pure fixed-node distributions in polynomial time, as well as not having
population fluctuations and biases. The current implementation does not
work with T-moves.

There are two adjustable parameters that affect the quality of the RMC
projection: imaginary projection time :math:`\beta` of the sampling path
(commonly called a “reptile") and the Trotter time step :math:`\tau`.
:math:`\beta` must be chosen to be large enough such that
:math:`e^{-\beta \hat{H}}|\Psi_T\rangle \approx |\Phi_0\rangle` for
mixed observables, and
:math:`e^{-\frac{\beta}{2} \hat{H}}|\Psi_T\rangle \approx |\Phi_0\rangle`
for pure observables. The reptile is discretized into
:math:`M=\beta/\tau` beads at the cost of an :math:`\mathcal{O}(\tau)`
time-step error for observables arising from the Trotter-Suzuki breakup
of the short-time propagator.

The following table lists some of the more practical

``vmc`` method:

  parameters:

  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | **Name**        | **Datatype** | **Values**              | **Default** | **Description**                                                |
  +=================+==============+=========================+=============+================================================================+
  | ``beta``        | real         | :math:`> 0`             | dep.        | Reptile project time :math:`\beta`                             |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``timestep``    | real         | :math:`> 0`             | 0.1         | Trotter time step :math:`\tau` for each electron move          |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``beads``       | int          | :math:`> 0`             | 1           | Number of reptile beads :math:`M=\beta/\tau`                   |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``blocks``      | integer      | :math:`> 0`             | 1           | Number of blocks                                               |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``steps``       | integer      | :math:`\geq 0`          | 1           | Number of steps per block                                      |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``vmcpresteps`` | integer      | :math:`\geq 0`          | 0           | Propagates reptile using VMC for given number of steps         |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``warmupsteps`` | integer      | :math:`\geq 0`          | 0           | Number of steps for warming up                                 |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+
  | ``maxAge``      | integer      | :math:`\geq 0`          | 0           | Force accept for stuck reptile if age exceeds ``maxAge``       |
  +-----------------+--------------+-------------------------+-------------+----------------------------------------------------------------+

Additional information:

Because of the sampling differences between DMC ensembles of walkers and
RMC reptiles, the RMC block should contain the following estimator
declaration to ensure correct sampling:
``<estimator name="RMC" hdf5="no">``.

-  ``beta`` or ``beads``? One or the other can be specified, and from
   the Trotter time step, the code will construct an appropriately sized
   reptile. If both are given, ``beta`` overrides ``beads``.

-  **Mixed vs. pure observables?** Configurations sampled by the
   endpoints of the reptile are distributed according to the mixed
   distribution
   :math:`f(\mathbf{R})=\Psi_T(\mathbf{R})\Phi_0(\mathbf{R})`. Any
   observable that is computable within DMC and is dumped to the
   ``scalar.dat`` file will likewise be found in the ``scalar.dat`` file
   generated by RMC, except there will be an appended ``_m`` to alert
   the user that the observable was computed on the mixed distribution.
   For pure observables, care must be taken in the interpretation. If
   the observable is diagonal in the position basis (in layman’s terms,
   if it is entirely computable from a single electron configuration
   :math:`\mathbf{R}`, like the potential energy), and if the observable
   does not have an explicit dependence on the trial wavefunction (e.g.,
   the local energy has an explicit dependence on the trial wavefunction
   from the kinetic energy term), then pure estimates will be correctly
   computed. These observables will be found in either the
   ``scalar.dat``, where they will be appended with a ``_p`` suffix, or
   in the ``stat.h5`` file. No mixed estimators will be dumped to the h5
   file.

-  **Sampling**: For pure estimators, the traces of both pure and mixed
   estimates should be checked. Ergodicity is a known problem in RMC.
   Because we use the bounce algorithm, it is possible for the reptile
   to bounce back and forth without changing the electron coordinates of
   the central beads. This might not easily show up with mixed
   estimators, since these are accumulated at constantly regrown ends,
   but pure estimates are accumulated on these central beads and so can
   exhibit strong autocorrelations in pure estimate traces.

-  **Propagator**: Our implementation of RMC uses Moroni’s DMC link
   action (symmetrized), with Umrigar’s scaled drift near nodes. In this
   regard, the propagator is identical to the one QMCPACK uses in DMC.

-  **Sampling**: We use Ceperley’s bounce algorithm. ``MaxAge`` is used
   in case the reptile gets stuck, at which point the code forces move
   acceptance, stops accumulating statistics, and requilibrates the
   reptile. Very rarely will this be required. For move proposals, we
   use particle-by-particle VMC a total of :math:`N_e` times to generate
   a new all-electron configuration, at which point the action is
   computed and the move is either accepted or rejected.

.. bibliography:: /bibs/methods.bib
