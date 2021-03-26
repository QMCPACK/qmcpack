#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

SET(test_utilities_src
  RandomForTest.cpp)
add_library(utilities_for_test ${test_utilities_src})
TARGET_INCLUDE_DIRECTORIES(utilities_for_test PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
