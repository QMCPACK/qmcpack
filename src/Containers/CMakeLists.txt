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

set(CONTAINER_SRCS MinimalContainers/ConstantSizeMatrix.cpp)
add_library(containers ${CONTAINER_SRCS})

target_include_directories(containers PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(containers PUBLIC platform_runtime)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
  add_subdirectory(OhmmsPETE/tests)
  add_subdirectory(OhmmsSoA/tests)
  add_subdirectory(MinimalContainers/tests)
  add_subdirectory(Pools/tests)
endif()
