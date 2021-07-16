#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(CPU_SRCS BlasThreadingEnv.cpp)

add_library(platform_cpu ${CPU_SRCS})
target_link_libraries(platform_cpu Math::BLAS_LAPACK Math::scalar_vector_functions)
