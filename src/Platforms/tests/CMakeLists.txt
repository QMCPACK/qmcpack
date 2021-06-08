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

if(BUILD_UNIT_TESTS)
  if(ENABLE_CUDA)
    add_subdirectory(CUDA)
  endif()

  add_subdirectory(OMPTarget)
  add_subdirectory(CPU)
endif()
