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

#ifndef AFQMC_DEVICE_KERNELS_HPP
#define AFQMC_DEVICE_KERNELS_HPP


#if defined(ENABLE_CUDA)

#include "AFQMC/Numerics/detail/CUDA/Kernels/determinant.cuh" 
#include "AFQMC/Numerics/detail/CUDA/Kernels/adotpby.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/fill_n.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/uninitialized_fill_n.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/uninitialized_copy_n.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/axty.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/adiagApy.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/sum.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/acAxpbB.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/print.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/setIdentity.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/zero_complex_part.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/batchedDot.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/copy_n_cast.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/inplace_cast.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/ajw_to_waj.cuh" 
#include "AFQMC/Numerics/detail/CUDA/Kernels/vKKwij_to_vwKiKj.cuh" 
#include "AFQMC/Numerics/detail/CUDA/Kernels/KaKjw_to_QKajw.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/vbias_from_v1.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/KaKjw_to_KKwaj.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/batched_dot_wabn_wban.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/batched_Tab_to_Klr.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/dot_wabn.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/Tab_to_Kl.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/sampleGaussianRNG.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/construct_X.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/reference_operations.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/term_by_term_matrix_vec.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/axpyBatched.cuh"

#endif


#endif
