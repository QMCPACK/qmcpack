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


IF (BUILD_UNIT_TESTS)
  IF(ENABLE_CUDA)
    SUBDIRS(CUDA)
  ENDIF()

  SUBDIRS(OMPTarget)
  SUBDIRS(CPU)
ENDIF()
