#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////
    
INCLUDE("${qmcpack_SOURCE_DIR}/CMake/macros.cmake")
    

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR wavefunction)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME unit_test_${SRC_DIR})

SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/wavefunction)
#SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT2 ${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT3 ${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/LiH-arb.pwscf.h5)

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT2} ${UTEST_DIR}/bccH.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT3} ${UTEST_DIR}/LiH-arb.pwscf.h5)

ADD_EXECUTABLE(${UTEST_EXE} test_wf.cpp test_bspline_jastrow.cpp test_einset.cpp test_pw.cpp
               test_polynomial_eeI_jastrow.cpp test_dirac_det.cpp test_dirac_matrix.cpp)
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmc qmcwfs qmcbase qmcutil ${QMC_UTIL_LIBS} ${MPI_LIBRARY})

ADD_UNIT_TEST(${UTEST_NAME} "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS "unit" WORKING_DIRECTORY ${UTEST_DIR})
