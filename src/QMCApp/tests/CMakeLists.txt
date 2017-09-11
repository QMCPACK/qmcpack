#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

# Also includes tests for code OhmmsApp directory

SET(SRC_DIR qmcapp)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})


SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/wavefunction)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.pwscf.h5)


ADD_EXECUTABLE(${UTEST_EXE} test_project_data.cpp test_rng_control.cpp test_particle_pool.cpp
                            test_wavefunction_pool.cpp)
#USE_FAKE_RNG(${UTEST_EXE})
#TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcdriver_unit qmcham qmcwfs qmcbase qmcutil qmcfakerng ${QMC_UTIL_LIBS} ${MPI_LIBRARY})
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcdriver_unit qmcham qmcwfs qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

#ADD_TEST(NAME ${UTEST_NAME} COMMAND "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit" WORKING_DIRECTORY ${UTEST_DIR})
