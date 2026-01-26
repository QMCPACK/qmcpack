#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

add_library(minimal_particle_pool MinimalParticlePool.cpp)
target_include_directories(minimal_particle_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_particle_pool PUBLIC qmcparticle)

set(SRC_DIR particle)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_${UTEST_EXE})

add_executable(
  ${UTEST_EXE}
  test_ParticleSet.cpp
  test_VirtualParticleSet.cpp
  test_walker.cpp
  test_particle_pool.cpp
  test_sample_stack.cpp
  test_DTModes.cpp
  test_MCCoords.cpp)
target_link_libraries(${UTEST_EXE} catch_main qmcparticle minimal_particle_pool)
if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} qmcutil qmcparticle_omptarget)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)

set(UTEST_EXE test_distance_table)
set(UTEST_NAME deterministic-unit_${UTEST_EXE})

add_executable(
  ${UTEST_EXE}
  test_distance_table.cpp
  test_SoaDistanceTableAA.cpp)
target_link_libraries(${UTEST_EXE} catch_main qmcparticle minimal_particle_pool)
if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} qmcutil qmcparticle_omptarget)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
