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

#ifndef AFQMC_ARCH_HPP
#define AFQMC_ARCH_HPP

#if defined(ENABLE_CUDA)
#include "AFQMC/Memory/CUDA/cuda_arch.h"
#elif defined(ENABLE_HIP)
#include "AFQMC/Memory/HIP/hip_arch.h"
#endif

#endif
