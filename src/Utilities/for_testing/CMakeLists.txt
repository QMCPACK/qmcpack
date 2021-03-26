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

if(USE_OBJECT_TARGET)
  add_library(utilities_for_test OBJECT ${test_utilities_src})
else()
  add_library(utilities_for_test ${test_utilties_src})
endif()

TARGET_INCLUDE_DIRECTORIES(utilities_for_test PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

  
