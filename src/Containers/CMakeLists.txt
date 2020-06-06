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


SET(CONTAINER_SRCS
  MinimalContainers/ConstantSizeMatrix.cpp
)
ADD_LIBRARY(containers ${CONTAINER_SRCS})

TARGET_INCLUDE_DIRECTORIES(containers PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_LINK_LIBRARIES(containers PUBLIC platform_device)

IF (BUILD_UNIT_TESTS)
  SUBDIRS(OhmmsPETE/tests)
  SUBDIRS(OhmmsSoA/tests)
  SUBDIRS(MinimalContainers/tests)
ENDIF()
