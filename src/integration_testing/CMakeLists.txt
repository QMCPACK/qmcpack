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
target_link_libraries(mock_walker_elements PUBLIC minimal_testing_pools qmcparticle qmcwfs qmcham message minimal_particle_pool)
target_link_libraries(mock_walker_elements PRIVATE catch_main)

  # a target for the minimal testing pools.
  # This should be used by test programs using any of the minimal
  # pools
  #
  # What it does:
  # It removes the need to include a relative
  # path into another modules testing directory to use a minimal pools
  #
  # What it doesn't do:
  # Provide linkage to the parent modules i.e. Particle, QMCHamiltonians,
  # QMCWaveFunction.
  #
  # It may be used faciliate transition to single compilation of the
  # minimal pools without a great deal of mess in the various test
  # program builds however it is useful as is. Consult with PD before
  # reverting this.
  #

add_library(minimal_testing_pools INTERFACE)
target_sources(minimal_testing_pools INTERFACE MinimalParticlePool.h MinimalHamiltonianPool.h MinimalWaveFunctionPool.h)
target_include_directories(minimal_testing_pools INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(minimal_particle_pool MinimalParticlePool.cpp)
target_include_directories(minimal_particle_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_particle_pool PUBLIC qmcparticle)
