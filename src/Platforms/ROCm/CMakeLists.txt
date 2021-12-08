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


add_library(platform_rocm_runtime INTERFACE)
target_link_libraries(platform_rocm_runtime INTERFACE hip::host)

set(ROCM_LA_SRCS hipBLAS.hip.cpp)
#set_source_files_properties(${ROCM_LA_SRCS} PROPERTIES HIP_SOURCE_PROPERTY_FORMAT 1)
set_source_files_properties(${ROCM_LA_SRCS} PROPERTIES LANGUAGE HIP)
add_library(platform_rocm_LA ${ROCM_LA_SRCS})
target_link_libraries(platform_rocm_LA PUBLIC roc::hipblas roc::rocsolver)
set_target_properties(platform_rocm_LA PROPERTIES LINKER_LANGUAGE HIP)
target_include_directories(platform_rocm_LA PUBLIC ${PROJECT_SOURCE_DIR}/external_codes/thrust)
target_link_libraries(platform_rocm_LA PRIVATE platform_rocm_runtime)
