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

set(WFSSRCS_CUDA determinant_update.cu delayed_update.cu PhaseFactors.cu AtomicOrbitalCuda.cu BsplineJastrowCuda.cu
                 BsplineJastrowCudaPBC.cu)

add_library(qmcwfs_cuda_legacy ${WFSSRCS_CUDA})
if(QMC_CUDA2HIP)
  set_source_files_properties(${WFSSRCS_CUDA} PROPERTIES LANGUAGE HIP)
endif()

target_link_libraries(qmcwfs_cuda_legacy PRIVATE platform_runtime qmc_external_thrust)
