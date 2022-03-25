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

The batched drivers introduce a new concept, "crowd", as a sub-organization of walker population.
A crowd is a subset of the walkers that are operated on as as single batch.
Walkers within a crowd operate their computation in lock-step, which helps the GPU efficiency.
Walkers between crowds remain fully asynchronous unless operations involving the full population are needed.
With this flexible batching capability the new drivers are capable of delivering maximal performance on given hardware.
In the new driver design, all the batched API calls may fallback to an existing single walker implementation.
Consequently, batched drivers allow mixing and matching CPU-only and GPU-accelerated features
in a way that is not feasible with the legacy GPU implementation.

For OpenMP GPU offload users, batched drivers are essential to effectively use GPUs.

