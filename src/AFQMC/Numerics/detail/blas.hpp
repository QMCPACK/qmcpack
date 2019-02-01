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

#ifndef AFQMC_BLAS_OPTIONS_HPP
#define AFQMC_BLAS_OPTIONS_HPP

#include<cassert>
#include "AFQMC/Numerics/detail/utilities.hpp"
#include "AFQMC/Numerics/detail/blas_cpu.hpp"
#if defined(QMC_CUDA)
#include "AFQMC/Numerics/detail/CUDA/blas_cuda.hpp"
#endif

#endif
