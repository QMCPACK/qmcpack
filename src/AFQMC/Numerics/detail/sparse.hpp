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

#ifndef AFQMC_SPARSE_OPTIONS_HPP
#define AFQMC_SPARSE_OPTIONS_HPP

#include <cassert>
#include "AFQMC/Numerics/detail/utilities.hpp"
#include "AFQMC/Numerics/detail/CPU/sparse_cpu.hpp"
#if defined(ENABLE_CUDA)
#include "AFQMC/Numerics/detail/CUDA/sparse_cuda_gpu_ptr.hpp"
#include "AFQMC/Numerics/detail/CUDA/sparse_cuda_catch_all.hpp"
#elif defined(ENABLE_HIP)
#include "AFQMC/Numerics/detail/HIP/sparse_hip_gpu_ptr.hpp"
#include "AFQMC/Numerics/detail/HIP/sparse_hip_catch_all.hpp"
#endif

#endif
