QMCPACK miniapps
================

Add `-DQMC_BUILD_LEVEL=5` in CMake to build the miniapps in this folder.
Optionally using `-DQMC_MIXED_PRECISION=1` to access mixed precision (most of the element values are single).
Assume that the mixed precision is sufficient to evaluate the algorithms, implementations and performance.


# checkpoint/restart miniapp
Its source and binary file are src/miniapps/restart.cpp and bin/restart.
It dumps walker configurations and random number seeds to the HDF5 files and then reads them in and check the correctness.
Pass or Fail will be printed at the end of the standard output.

Parallel Collective I/O is implemented via parallel HDF5. It is enabled by default when parallel HDF5 library is available.
To have good performance at large scale, version 1.10 is needed.
