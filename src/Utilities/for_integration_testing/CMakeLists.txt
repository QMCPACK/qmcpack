#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

add_library(integration_testing MockGoldWalkerElements.cpp)
target_include_directories(integration_testing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(integration_testing PUBLIC minimal_testing_pools qmcparticle qmcwfs message)
target_link_libraries(integration_testing PRIVATE catch_main)
