#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

# a targets for the individual minimal testing pools.
#
# They may be used to minimize dependency in tests when possible.
#

add_library(minimal_particle_pool MinimalParticlePool.cpp)
target_include_directories(minimal_particle_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_particle_pool PUBLIC qmcparticle)

add_library(minimal_wavefunction_pool MinimalWaveFunctionPool.cpp)
target_include_directories(minimal_wavefunction_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_wavefunction_pool PUBLIC qmcwfs)

add_library(minimal_hamiltonian_pool MinimalHamiltonianPool.cpp)
target_include_directories(minimal_hamiltonian_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_hamiltonian_pool PUBLIC qmcham)

# use this target if using all three minimal pools
add_library(minimal_test_pools INTERFACE)
target_link_libraries(minimal_test_pools INTERFACE minimal_particle_pool minimal_wavefunction_pool minimal_hamiltonian_pool)

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
target_link_libraries(mock_walker_elements PUBLIC qmcparticle qmcwfs qmcham minimal_test_pools message)
target_link_libraries(mock_walker_elements PRIVATE catch_main)
