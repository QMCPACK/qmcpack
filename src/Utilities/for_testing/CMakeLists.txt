#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

SET(test_utilities_src
  RandomForTest.cpp
  checkMatrix.cpp)
add_library(utilities_for_test ${test_utilities_src})
target_include_directories(utilities_for_test PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(utilities_for_test PRIVATE "${PROJECT_SOURCE_DIR}/external_codes/catch")

add_library(catch_main_no_mpi catch_main_no_mpi.cpp)
TARGET_INCLUDE_DIRECTORIES(catch_main_no_mpi PUBLIC "${PROJECT_SOURCE_DIR}/external_codes/catch")
