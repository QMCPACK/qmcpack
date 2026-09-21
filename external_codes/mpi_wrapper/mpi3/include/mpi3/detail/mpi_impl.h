// Copyright 2025-2026 Alfredo A. Correa

#define OMPI_SKIP_MPICXX 1  // NOLINT  // https://github.com/open-mpi/ompi/issues/5157
// #if __has_include(<mpi.h>)
#include<mpi.h>
// #else
// #pragma message("since <mpi.h> was not found, the library is falling back to an internal version of the header <mpi3/mpich_hack/mpi.h>, this is to check for compilation but your build may not link (or execute incorrect code)")
// #include <mpi3/mpich_hack/mpi.h>
// #endif
