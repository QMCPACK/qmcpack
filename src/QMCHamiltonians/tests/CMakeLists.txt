#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////
    
    

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR hamiltonian)
SET(UTEST_EXE test_${SRC_DIR})
SET(UTEST_NAME deterministic-unit_test_${SRC_DIR})
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/hamiltonians)
SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT} ${UTEST_DIR}/pwscf.pwscf.h5)


SET(SRCS test_bare_kinetic.cpp
         test_coulomb_pbcAB.cpp
         test_coulomb_pbcAB_ewald.cpp
         test_coulomb_pbcAA.cpp
         test_coulomb_pbcAA_ewald.cpp
         test_force.cpp
         test_force_ewald.cpp
         test_stress.cpp
         test_spacewarp.cpp
         test_ecp.cpp
         test_hamiltonian_pool.cpp
         test_hamiltonian_factory.cpp
         test_PairCorrEstimator.cpp
         test_SkAllEstimator.cpp
         test_QMCHamiltonian.cpp
         test_derivs.cpp
         test_ObservableHelper.cpp
         )
         
IF(QMC_CUDA)
  SET(SRCS ${SRCS}
      test_coulomb_CUDA.cpp 
)
ENDIF()

EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")

foreach(fname Na2.structure.xml simple.txt)
  file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/${fname} ${UTEST_DIR}/${fname} SYMBOLIC)
endforeach()

foreach(fname cn.wfnoj.xml cn.wfj.xml cn.msd-wfnoj.xml cn.msd-wfj.xml)
  file(CREATE_LINK ${qmcpack_SOURCE_DIR}/src/QMCWaveFunctions/tests/${fname} ${UTEST_DIR}/${fname} SYMBOLIC)
endforeach()

foreach(fname C.BFD.xml Na.BFD.xml so_ecp_test.xml C.ccECP.xml N.ccECP.xml)
  file(CREATE_LINK ${qmcpack_SOURCE_DIR}/tests/pseudopotentials_for_tests/${fname} ${UTEST_DIR}/${fname} SYMBOLIC)
endforeach()

ADD_EXECUTABLE(${UTEST_EXE} ${SRCS})
TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcham)
IF(USE_OBJECT_TARGET)
TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcwfs qmcparticle qmcutil containers platform_omptarget)
ENDIF()

ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

IF( "${ENABLE_SANITIZER}" STREQUAL "asan" ) 
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS noasan)
ENDIF()
