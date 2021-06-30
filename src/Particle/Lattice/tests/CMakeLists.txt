#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

set(SRC_DIR lattice)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})

add_executable(${UTEST_EXE} test_ParticleBConds.cpp test_CrystalLattice.cpp test_LRBreakupParameters.cpp)
target_include_directories(${UTEST_EXE} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../..")
target_link_libraries(${UTEST_EXE} catch_main qmcutil Math::scalar_vector_functions)

add_unit_test(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
