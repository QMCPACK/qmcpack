#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(ROCM_SRCS hipBLAS.hip.cpp)

set_source_files_properties(${ROCM_SRCS} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
hip_add_library(platform_rocm ${ROCM_SRCS})
target_link_libraries(platform_rocm PUBLIC hip::host roc::hipblas roc::rocsolver)
target_include_directories(platform_rocm PUBLIC ${PROJECT_SOURCE_DIR}/external_codes/thrust)
