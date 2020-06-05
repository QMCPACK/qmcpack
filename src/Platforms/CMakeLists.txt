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
    sysutil.cpp
    Host/InfoStream.cpp
    Host/OutputManager.cpp
)

SET(DEVICE_SRCS
    accelerators.cpp
)

ADD_LIBRARY(platform_host ${HOST_SRCS})

ADD_LIBRARY(platform_device ${DEVICE_SRCS})
TARGET_INCLUDE_DIRECTORIES(platform_device PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_LINK_LIBRARIES(platform_device PUBLIC message)

SUBDIRS(OpenMP)
TARGET_LINK_LIBRARIES(platform_device PRIVATE platform_omp)

IF(ENABLE_CUDA)
  SUBDIRS(CUDA)
  TARGET_LINK_LIBRARIES(platform_device PUBLIC platform_cuda)
ENDIF(ENABLE_CUDA)
