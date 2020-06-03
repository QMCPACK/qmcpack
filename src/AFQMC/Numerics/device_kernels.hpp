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
#include "AFQMC/Numerics/detail/CUDA/Kernels/Auwn_Bun_Cuw.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/inplace_product.cuh"
#include "AFQMC/Numerics/detail/CUDA/Kernels/get_diagonal.cuh"

#elif defined(ENABLE_HIP)

#include "AFQMC/Numerics/detail/HIP/Kernels/determinant.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/adotpby.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/fill_n.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/uninitialized_fill_n.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/uninitialized_copy_n.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/axty.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/adiagApy.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/sum.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/acAxpbB.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/print.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/setIdentity.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/zero_complex_part.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/batchedDot.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/copy_n_cast.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/inplace_cast.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/ajw_to_waj.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/vKKwij_to_vwKiKj.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/KaKjw_to_QKajw.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/vbias_from_v1.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/KaKjw_to_KKwaj.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/batched_dot_wabn_wban.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/batched_Tab_to_Klr.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/dot_wabn.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/Tab_to_Kl.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/sampleGaussianRNG.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/construct_X.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/reference_operations.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/term_by_term_matrix_vec.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/axpyBatched.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/Auwn_Bun_Cuw.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/inplace_product.hip.h"
#include "AFQMC/Numerics/detail/HIP/Kernels/get_diagonal.hip.h"

#endif

#endif
