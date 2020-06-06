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


SET(CPU_SRCS
)

ADD_LIBRARY(platform_cpu INTERFACE)
TARGET_INCLUDE_DIRECTORIES(platform_cpu INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

IF (BUILD_UNIT_TESTS)
  SUBDIRS(tests)
ENDIF()
