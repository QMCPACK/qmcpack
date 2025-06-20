#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

# at the crowd thread scope of the API's each walker is associated
# with one object of ParticleSet, QMCWavefunction, and
# QMCHamiltonian i.e. walker elements
# This integration library provides helpers for integration tests
# working at this API level
#
# Right now it only provides a helper class and factory functions for
# the gold elements from which the individual walker elements are
# cloned but I expect refactor repeated code for the vectors of walker
# elements needed for many integration tests in the future
add_library(mock_walker_elements MockGoldWalkerElements.cpp)
target_include_directories(mock_walker_elements PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(mock_walker_elements PUBLIC minimal_testing_pools qmcparticle qmcwfs qmcham message)
target_link_libraries(mock_walker_elements PRIVATE catch_main)
