.. _performance_portable:

Performance Portable Implementation
===================================

Under the Exascale Computing Project effort a new set of QMC drivers was developed
to eliminate the divergence of legacy CPU and GPU code paths at the QMC driver level and make the drivers CPU/GPU agnostic.
The divergence came from the the fact that the CPU code path favors executing all the compute tasks within a step
for one walker and then advance walker by walker. Multiple CPU threads process their own assigned walkers in parallel.
In this way, walkers are not synchronized with each other and maximal throughout can be achieved on CPU.
The GPU code path favors executing the same compute task over all the walkers together to maximize GPU throughput.
This compute dispatch pattern minimizes the overhead of dispatching computation and host-device data transfer.
However, the legacy GPU code path only leverages the OpenMP main host thread for handling
all the interaction between the host and GPUs and limit the kernel dispatch capability.
In brief, the CPU code path handles computation with a walker batch size of one and many batches
while the GPU code path uses only one batch containing all the walkers.
The new drivers that implement this flexible batching scheme are called "batched drivers".

For OpenMP GPU offload users, batched drivers are essential to effectively use GPUs.

Links to more information in other sections of the manual:

 - **Build instructions:** :ref:`OpenMP target offload <offloadbuild>` section of the :ref:`obtaininginstalling` chapter.

 - **Supported features:** :ref:`gpufeatures` section of the :ref:`chap:features` chapter.

 - **Driver Inputs:** :ref:`batched_drivers` section of the :ref:`qmcmethods` chapter.
