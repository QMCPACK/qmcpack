#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


# target for host and programming model runtime systems
set(DEVICE_SRCS MemoryUsage.cpp DualAllocator.cpp DeviceManager.cpp)
add_library(platform_runtime ${DEVICE_SRCS})
target_link_libraries(platform_runtime PUBLIC platform_host_runtime)
target_include_directories(platform_runtime PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# target for host and device linear algebra libraries
add_library(platform_LA INTERFACE)
target_link_libraries(platform_LA INTERFACE platform_runtime)

set(HOST_SRCS Host/sysutil.cpp Host/InfoStream.cpp Host/OutputManager.cpp)
add_library(platform_host_runtime ${HOST_SRCS})

add_subdirectory(CPU)
target_link_libraries(platform_runtime PUBLIC platform_cpu_runtime)
target_link_libraries(platform_LA INTERFACE platform_cpu_LA)

add_subdirectory(OMPTarget)
target_link_libraries(platform_runtime PUBLIC platform_omptarget_runtime)
target_link_libraries(platform_LA INTERFACE platform_omptarget_LA)

if(QMC_CUDA)
  add_subdirectory(CUDA_legacy)
  target_link_libraries(platform_runtime PUBLIC platform_cuda_legacy)
endif(QMC_CUDA)

if(ENABLE_CUDA)
  add_subdirectory(CUDA)
  target_link_libraries(platform_runtime PUBLIC platform_cuda_runtime)
  target_link_libraries(platform_LA INTERFACE platform_cuda_LA)
endif(ENABLE_CUDA)

if(ENABLE_ROCM)
  add_subdirectory(ROCm)
  target_link_libraries(platform_runtime PUBLIC platform_rocm_runtime)
  target_link_libraries(platform_LA INTERFACE platform_rocm_LA)
endif(ENABLE_ROCM)

add_subdirectory(tests)
