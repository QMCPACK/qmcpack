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

set(NUMERICS_SRCS SmoothFunctions.cpp OneDimGridFactory.cpp SoaSphericalTensor.cpp)

add_library(qmcnumerics ${NUMERICS_SRCS})
target_link_libraries(qmcnumerics qmcio containers)
target_include_directories(qmcnumerics INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
