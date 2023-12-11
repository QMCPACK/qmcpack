.. _performance_portable:

Performance Portable Implementation
===================================

The so-called performance portable implementation was developed to present a unified way to run QMC on CPU and GPU
systems, and eliminate the divergence between CPU and GPU code paths that had been introduced in the past, while still
maintaining high performance. This required generalizing all the driver inputs to potentially drive larger batches of
walkers and also eliminating ambiguities in the various input blocks of QMCPACK. Internally many new code paths have
been created, including new QMC drivers for VMC, DMC, and the wavefunction optimizer. 

Once this implementation is sufficiently matured and enough features are available, the old non-performance portable
drivers will be deprecated and eventually deleted. The number of changes required to old input files is usually very
small, so use of the new performance portable implementation is encouraged, particularly for new projects.

The performance portable implementation load balances the total number of walkers onto MPI tasks, as per the old
drivers. The new implementation is then able to subdivide the walkers of each MPI task into multiple similarly-sized
crowds. The walkers in each crowd can then be updated simultaneously. This structure enables the walkers to be
efficiently mapped to both CPUs and GPUs. On CPU systems, they then are mapped to OpenMP threads where a single walker
can be computed efficiently by even a single thread. On GPU systems, large numbers of GPU threads must be used
concurrently for high efficiency: Each crowd is first owned by a distinct CPU thread, which in turn executes batched
operations over all the walkers in its crowd on the GPU. Provided the batches are sufficiently large, this facilitates
efficient GPU execution, while the use of multiple crowds can reduce synchronization and allow higher performance to be
obtained. For these reasons the new performance portable drivers are also referred to as batched drivers, since this is
the largest change from the older code.

The new implementation largely uses OpenMP offload for portability, although other technologies are also used and the
implementation has flexible dispatch to help obtain high performance on every platform.

This implementation was designed and implemented as part of the Exascale Computing Project, with a view to bringing
QMCPACK to GPUs from multiple vendors with high-efficiency while creating a more maintainable and easy to contribute to
codebase.

Links to more information in other sections of the manual:

 - **Build instructions:** :ref:`OpenMP target offload <offloadbuild>` section of the :ref:`obtaininginstalling` chapter.

 - **Supported features:** :ref:`gpufeatures` section of the :ref:`chap:features` chapter.

 - **Enabling batch drivers** :ref:`driver-version-parameter` section of the :ref:`input-overview` chapter.

 - **Driver Inputs:** :ref:`batched_drivers` section of the :ref:`qmcmethods` chapter.


Input files for batched drivers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following changes to update input files to use the batched drivers.

1. Update the project block with the ``driver_version`` parameter. For example:

::

  <project id="vmc" series="0">
     <parameter name="driver_version">batch</parameter>
  </project>

See :ref:`driver-version-parameter` for more.

2. Modify the QMC algorithm blocks

The most significant change is the ``walkers`` parameter has been replaced with ``walkers_per_rank`` or ``total_walkers``.

See  :ref:`batched_drivers` for details.
