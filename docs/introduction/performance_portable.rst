.. _performance_portable:

QMCPACK's Performance Portable Implementation
=============================================

With the release of QMCPACK v4.0, the default execution paths utilize the so-called performance portable design and implementation,
also sometimes referred to as the batched code or batched drivers. The implementation includes new drivers for each QMC algorithm through
to updated wavefunctions and observables. This implementation was developed to present a unified way to run QMC on CPU and GPU
systems, and eliminate the divergence between CPU and GPU code paths that had been introduced in the past. In practice this was hard to
maintain and the feature sets were never equal. Performance is maintained or increased in the new design. 
Note that these updates required generalizing all the driver inputs to drive batches of walkers and eliminating historical ambiguities
in the various input blocks of QMCPACK. Consequently, small updates to old input files may be required.

The performance portable implementation is now considered sufficiently mature and feature rich to support the majority of
calculations done with QMCPACK. There is no CMake configuration required to enable this implementation -- all builds include it.
The legacy drivers and implementations are considered deprecated and will be removed
in upcoming QMCPACK versions. Their removal and the subsequent reduction in behind-the-scenes entangled data structures and duplication of
similar code will enable the implementation of more efficient QMC algorithms and allow the development team to maximize their productivity.

We do not recommend starting any new projects using the old non-performance portable legacy drivers or code.

The performance portable implementation load balances the total number of walkers onto MPI tasks. While the old DMC driver could do 
this, in the new implementation even the VMC driver is able to drive large numbers of walkers.
The new implementation is then able to subdivide the walkers of each MPI task into multiple similarly-sized
crowds. This is a new level of parallelization that can enable higher performance. The walkers in each crowd can then be updated
simultaneously. This structure enables the walkers to be efficiently mapped to both CPUs and GPUs. On CPU systems, they then are mapped to OpenMP threads
where individual walker can be updated efficiently by a single thread. On GPU systems, large numbers of GPU threads must be used
concurrently for high efficiency: Each crowd is first owned by a distinct CPU thread, which in turn executes batched
operations over all the walkers in its crowd on the GPU. Provided the batches are sufficiently large, this facilitates
efficient GPU execution, while the use of multiple crowds can reduce synchronization and allow higher performance to be
obtained. For these reasons the new performance portable drivers are often referred to as batched drivers, since this is
the largest change from the older code.

The new implementation largely uses OpenMP offload for portability, although other technologies are also used and the
implementation has flexible dispatch to help obtain high performance on every platform. If needed, custom highly optimized implementations
could be added for specific architectures or architectural features without a substantial rewrite.

QMCPACK currently fully GPU accelerates spline and LCAO orbitals, single and multideterminant wavefunctions, all-electron and pseudopotential
calculations, one and two body Jastrow functions, open and closed boundary conditions, and the standard variational, wavefunction optimization,
and diffusion Monte Carlo algorithms. Other options and methods run via the CPU host. As a result, nearly all research-level molecular and
materials calculations achieve significant acceleration.

This implementation was designed and initially implemented as part of the Exascale Computing Project, with a view to bringing
QMCPACK to GPUs from multiple vendors with high-efficiency while creating a more maintainable and easy to contribute to
codebase.

Links to more information in other sections of the manual:

 - **GPU related build options:** :ref:`cmakeoptions` section of the :ref:`obtaininginstalling` chapter.

 - **Selecting between batch and legacy drivers** :ref:`driver-version-parameter` section of the :ref:`input-overview` chapter.

 - **Driver Inputs:** :ref:`batched_drivers` section of the :ref:`qmcmethods` chapter.


Updating input files for batched drivers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following changes to update input files to use the batched drivers.

1. Update the project block to include the ``driver_version`` parameter. While not required (their use is the default), specifying the ``driver_version``
   helps indicate that the input is for the modern code. For example:

::

  <project id="my_qmc_research" series="0">
     <parameter name="driver_version">batch</parameter>
  </project>

See :ref:`driver-version-parameter` for details.

2. Modify the QMC algorithm blocks

The most significant change is the ``walkers`` parameter has been replaced with ``walkers_per_rank`` or ``total_walkers``.

See :ref:`batched_drivers` for details.
