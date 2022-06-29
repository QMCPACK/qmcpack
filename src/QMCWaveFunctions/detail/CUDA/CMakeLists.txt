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

add_library(qmcwfs_cuda ${WFSSRCS_CUDA})
if(QMC_CUDA2HIP)
  set_source_files_properties(${WFSSRCS_CUDA} PROPERTIES LANGUAGE HIP)
endif()

target_link_libraries(qmcwfs_cuda PRIVATE platform_LA qmc_external_thrust)
