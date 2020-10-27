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

#ifndef AFQMC_CUSTOM_POINTERS_HPP
#define AFQMC_CUSTOM_POINTERS_HPP

#include "AFQMC/Memory/raw_pointers.hpp"
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
#include "AFQMC/Memory/device_pointers.hpp"
//#include "AFQMC/Memory/CUDA/cuda_utilities.h"
//#include "AFQMC/Memory/CUDA/cuda_init.h"
//#include "AFQMC/Memory/CUDA/cuda_gpu_pointer.hpp"
#endif

#endif
