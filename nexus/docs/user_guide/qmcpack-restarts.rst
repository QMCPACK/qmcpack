Running QMCPACK restarts
========================

Request a restart through a simulation dependency:

.. code-block:: python

  qmc2 = generate_qmcpack(
      ...,
      dependencies = [(conv,'orbitals'),
                      (qmc1,'restart')],
      )

The final QMC section in ``qmc1`` must enable checkpointing.  Nexus reads its
continuation XML and checks for both the walker and random-state HDF5 files.

Separate directories
--------------------

Give the two simulations different paths.  Nexus writes an ``mcwalkerset``
with a relative path to the checkpoint and retains the project ID of ``qmc2``.
See ``05_diamond_dft_dmc_restart/diamond_lda_dmc_restart_separate_dirs.py``
in ``examples/qmcpack/rsqmc_quantum_espresso`` for a complete workflow.

.. code-block:: python

  qmc1 = generate_qmcpack(identifier='dmc',path='diamond/dmc_initial',...)
  qmc2 = generate_qmcpack(
      identifier   = 'dmc_restart',
      path         = 'diamond/dmc_restart',
      dependencies = (qmc1,'restart'),
      ...,
      )

Same directory
--------------

Use distinct Nexus identifiers but the same path.  Nexus continues the
QMCPACK project ID and series from the continuation XML, producing successive
files such as ``dmc.s001.*`` and ``dmc.s002.*`` without colliding Nexus input,
output, error, or analyzer files.
See ``05_diamond_dft_dmc_restart/diamond_lda_dmc_restart_same_dir.py`` for a
complete workflow.

.. code-block:: python

  qmc1 = generate_qmcpack(identifier='dmc',path='diamond/dmc',...)
  qmc2 = generate_qmcpack(
      identifier   = 'dmc_restart',
      path         = 'diamond/dmc',
      dependencies = (qmc1,'restart'),
      ...,
      )

Twist averaging
---------------

The same patterns apply to twist-averaged calculations.  The initial and
restarted simulations must have matching twist grids; Nexus pairs and checks
the checkpoint for each twist.  Both directory layouts are demonstrated in
``examples/qmcpack/rsqmc_quantum_espresso/06_diamond_dft_dmc_twistavg_restart``.
