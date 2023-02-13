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

# Targets under Platforms are designed to abstract external libraries
# platform_runtime, platform_LA may be used outside Platforms.
# Device specific targets should not be referenced used outside Platforms.

# platform_runtime is for host and programming model runtime systems which inclues
# Device management: device assignement, memory management. Note: CPU is a device
# Math functions: scalar and vector math functions from OS or vendors
set(DEVICE_SRCS MemoryUsage.cpp DualAllocator.cpp DeviceManager.cpp PlatformSelector.cpp)
add_library(platform_runtime ${DEVICE_SRCS})
target_link_libraries(platform_runtime PUBLIC platform_host_runtime)
target_include_directories(platform_runtime PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# platform_LA is for host and device linear algebra libraries
add_library(platform_LA INTERFACE)
target_link_libraries(platform_LA INTERFACE platform_runtime)

# platform_host_runtime is the target for host runtime system which includes
# interaction with OS libraries
add_subdirectory(Host)

# include CPU platform
add_subdirectory(CPU)
target_link_libraries(platform_runtime PUBLIC platform_cpu_runtime)
target_link_libraries(platform_LA INTERFACE platform_cpu_LA)

# include OpenMP targert platform
add_subdirectory(OMPTarget)
target_link_libraries(platform_runtime PUBLIC platform_omptarget_runtime)
target_link_libraries(platform_LA INTERFACE platform_omptarget_LA)

# include CUDA platform
if(ENABLE_CUDA)
  add_subdirectory(CUDA)
  target_link_libraries(platform_runtime PUBLIC platform_cuda_runtime)
  target_link_libraries(platform_LA INTERFACE platform_cuda_LA)
endif(ENABLE_CUDA)

# include ROCm platform
if(ENABLE_ROCM)
  add_subdirectory(ROCm)
  target_link_libraries(platform_runtime PUBLIC platform_rocm_runtime)
  target_link_libraries(platform_LA INTERFACE platform_rocm_LA)
endif(ENABLE_ROCM)

# include SYCL platform
if(ENABLE_SYCL)
  add_subdirectory(SYCL)
  target_link_libraries(platform_runtime PUBLIC platform_sycl_runtime)
  target_link_libraries(platform_LA INTERFACE platform_sycl_LA)
endif(ENABLE_SYCL)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
