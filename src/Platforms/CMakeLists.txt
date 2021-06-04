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


SET(HOST_SRCS
    Host/sysutil.cpp
    Host/InfoStream.cpp
    Host/OutputManager.cpp
)

SET(DEVICE_SRCS
  MemoryUsage.cpp
  DualAllocator.cpp
)

ADD_LIBRARY(platform_host ${HOST_SRCS})

IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(platform_device OBJECT ${DEVICE_SRCS})
ELSE(USE_OBJECT_TARGET)
  ADD_LIBRARY(platform_device ${DEVICE_SRCS})
ENDIF(USE_OBJECT_TARGET)

TARGET_INCLUDE_DIRECTORIES(platform_device PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

SUBDIRS(CPU)
TARGET_LINK_LIBRARIES(platform_device PUBLIC platform_cpu)

SUBDIRS(OMPTarget)
TARGET_LINK_LIBRARIES(platform_device PUBLIC platform_omptarget)

IF(QMC_CUDA)
  SUBDIRS(CUDA_legacy)
  TARGET_LINK_LIBRARIES(platform_device PUBLIC platform_cuda_legacy)
ENDIF(QMC_CUDA)

IF(ENABLE_CUDA)
  SUBDIRS(CUDA)
  TARGET_LINK_LIBRARIES(platform_device PUBLIC platform_cuda)
ENDIF(ENABLE_CUDA)

SUBDIRS(tests)
