//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_LAPACK_OPTIONS_HPP
#define AFQMC_LAPACK_OPTIONS_HPP

#include <cassert>
#include "AFQMC/Numerics/detail/CPU/lapack_cpu.hpp"
#if defined(ENABLE_CUDA)
#include "AFQMC/Numerics/detail/CUDA/lapack_cuda_gpu_ptr.hpp"
//#include "AFQMC/Numerics/detail/CUDA/lapack_cuda_catch_all.hpp"
#elif defined(ENABLE_HIP)
#include "AFQMC/Numerics/detail/HIP/lapack_hip_gpu_ptr.hpp"
#endif
#endif
