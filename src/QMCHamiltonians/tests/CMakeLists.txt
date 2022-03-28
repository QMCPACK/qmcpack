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

set(SRC_DIR hamiltonian)
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(COULOMB_SRCS test_coulomb_pbcAB.cpp test_coulomb_pbcAB_ewald.cpp test_coulomb_pbcAA.cpp
                 test_coulomb_pbcAA_ewald.cpp test_EwaldRef.cpp)
set(FORCE_SRCS test_force.cpp test_force_ewald.cpp test_stress.cpp test_spacewarp.cpp)
set(HAM_SRCS
    test_bare_kinetic.cpp
    test_NonLocalTOperator.cpp
    test_ecp.cpp
    test_hamiltonian_pool.cpp
    test_hamiltonian_factory.cpp
    test_PairCorrEstimator.cpp
    test_SkAllEstimator.cpp
    test_QMCHamiltonian.cpp
    test_ObservableHelper.cpp)

if(QMC_CUDA)
  set(COULOMB_SRCS ${COULOMB_SRCS} test_coulomb_CUDA.cpp)
else()
  set(FORCE_SRCS ${FORCE_SRCS} test_ion_derivs.cpp)
endif()

set(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)

foreach(fname Na2.structure.xml simple.txt)
  maybe_symlink(${CMAKE_CURRENT_SOURCE_DIR}/${fname} ${UTEST_DIR}/${fname})
endforeach()

foreach(fname cn.wfnoj.xml cn.wfj.xml cn.msd-wfnoj.xml cn.msd-wfj.xml)
  maybe_symlink(${qmcpack_SOURCE_DIR}/src/QMCWaveFunctions/tests/${fname} ${UTEST_DIR}/${fname})
endforeach()

foreach(fname C.BFD.xml Na.BFD.xml so_ecp_test.xml C.ccECP.xml N.ccECP.xml)
  maybe_symlink(${qmcpack_SOURCE_DIR}/tests/pseudopotentials_for_tests/${fname} ${UTEST_DIR}/${fname})
endforeach()

foreach(CATEGORY coulomb force ham)
  set(UTEST_EXE test_${SRC_DIR}_${CATEGORY})
  set(UTEST_NAME deterministic-unit_${UTEST_EXE})
  string(TOUPPER "${CATEGORY}_SRCS" SOURCE_FILE_VAR_NAME)
  add_executable(${UTEST_EXE} ${${SOURCE_FILE_VAR_NAME}})

  target_link_libraries(${UTEST_EXE} catch_main qmcham)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} qmcwfs qmcparticle qmcutil platform_omptarget_LA)
  endif()

  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
endforeach()
