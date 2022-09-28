#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2022 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(WFSSRCS_CUDA cuBLAS_LU.cu delayed_update_helper.cu matrix_update_helper.cu)

if(NOT QMC_CUDA2HIP)
  add_library(qmcwfs_cuda ${WFSSRCS_CUDA})
else()
  hip_add_library(qmcwfs_cuda ${WFSSRCS_CUDA})
endif(NOT QMC_CUDA2HIP)

target_link_libraries(qmcwfs_cuda PRIVATE platform_LA qmc_external_thrust)
