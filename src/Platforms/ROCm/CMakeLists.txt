#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

add_library(platform_rocm_runtime INTERFACE)
target_link_libraries(platform_rocm_runtime INTERFACE hip::host)

set(ROCM_LA_SRCS hipBLAS.cpp)
add_library(platform_rocm_LA ${ROCM_LA_SRCS})
target_link_libraries(platform_rocm_LA PUBLIC roc::hipblas roc::rocsolver)
target_link_libraries(platform_rocm_LA PRIVATE platform_rocm_runtime)

# Add HIPBLAS_V2 for hipBLAS < 3.0.0
if(hipblas_VERSION VERSION_LESS "3.0.0")
  target_compile_definitions(platform_rocm_LA PUBLIC -DHIPBLAS_V2)
  message(STATUS "Defining HIPBLAS_V2 for compatibility with hipBLAS < 3.0.0 (Detected ${hipblas_VERSION})")
endif()
