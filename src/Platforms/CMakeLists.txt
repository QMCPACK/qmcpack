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


# Too much burden has been put on platform_device. It will be deprecated.
# We will move to functionality based targets
# platform_runtime handles device management, memory management
# add_library(platform_runtime INTERFACE)
# platform_lapack handles BLAS/LAPACK
# add_library(platform_lapack INTERFACE)

set(HOST_SRCS Host/sysutil.cpp Host/InfoStream.cpp Host/OutputManager.cpp)

set(DEVICE_SRCS MemoryUsage.cpp DualAllocator.cpp DeviceManager.cpp)

add_library(platform_host ${HOST_SRCS})

if(USE_OBJECT_TARGET)
  add_library(platform_device OBJECT ${DEVICE_SRCS})
else(USE_OBJECT_TARGET)
  add_library(platform_device ${DEVICE_SRCS})
endif(USE_OBJECT_TARGET)

target_include_directories(platform_device PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(platform_device PUBLIC platform_host)

add_subdirectory(CPU)
target_link_libraries(platform_device PUBLIC platform_cpu)

add_subdirectory(OMPTarget)
target_link_libraries(platform_device PUBLIC platform_omptarget)

if(QMC_CUDA)
  add_subdirectory(CUDA_legacy)
  target_link_libraries(platform_device PUBLIC platform_cuda_legacy)
endif(QMC_CUDA)

if(ENABLE_CUDA)
  add_subdirectory(CUDA)
  target_link_libraries(platform_device PUBLIC platform_cuda)
endif(ENABLE_CUDA)

add_subdirectory(tests)
