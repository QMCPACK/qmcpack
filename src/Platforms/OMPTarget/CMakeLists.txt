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


SET(OMP_SRCS
  ompBLAS.cpp
  OMPallocator.cpp
)

IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(platform_omptarget OBJECT ${OMP_SRCS})
ELSE()
  ADD_LIBRARY(platform_omptarget ${OMP_SRCS})
ENDIF()

IF(NOT QMC_COMPLEX)
  TARGET_COMPILE_DEFINITIONS(platform_omptarget PRIVATE -DOPENMP_NO_COMPLEX)
ENDIF()
