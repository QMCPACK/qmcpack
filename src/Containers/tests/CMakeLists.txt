#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(container_testing_src makeRngSpdMatrix.cpp)
add_library(container_testing ${container_testing_src})
target_include_directories(container_testing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(container_testing utilities_for_test containers)
