#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2025 QMCPACK developers.
#//
#// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
#//
#// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

add_library(minimal_hamiltonian_pool MinimalHamiltonianPool.cpp)
target_include_directories(minimal_hamiltonian_pool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(minimal_hamiltonian_pool PUBLIC qmcham)

set(SRC_DIR hamiltonian)
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(COULOMB_SRCS test_coulomb_pbcAB.cpp test_coulomb_pbcAB_ewald.cpp test_coulomb_pbcAA.cpp
                 test_coulomb_pbcAA_ewald.cpp test_EwaldRef.cpp test_NonLocalECPotential.cpp)
set(EWALD2D_SRCS test_ewald2d.cpp test_ewald_quasi2d.cpp)
set(FORCE_SRCS test_force.cpp test_force_ewald.cpp test_stress.cpp test_spacewarp.cpp)
set(HAM_SRCS
    test_bare_kinetic.cpp
    test_density_estimator.cpp
    test_L2Potential.cpp
    test_NonLocalTOperator.cpp
    test_ecp_builder.cpp
    test_ecp.cpp
    test_hamiltonian_pool.cpp
    test_hamiltonian_factory.cpp
    test_PairCorrEstimator.cpp
    test_SkAllEstimator.cpp
    test_QMCHamiltonian.cpp
    test_ObservableHelper.cpp
    test_Listener.cpp)

if(NOT QMC_COMPLEX)
  set(HAM_SRCS ${HAM_SRCS} test_RotatedSPOs_NLPP.cpp)
else()
  set(HAM_SRCS ${HAM_SRCS} test_SOECPotential.cpp)
endif()

set(FORCE_SRCS ${FORCE_SRCS} test_ion_derivs.cpp)

set(ALIASED_FILES)
make_file_alias(${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5 diamondC_1x1x1.pwscf.h5)
make_file_alias(${qmcpack_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp/pwscf.pwscf.h5 hcpBe.pwscf.h5)

set(UTEST_INPUT_FILES ${ALIASED_FILES}
    ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/C_diamond-twist-third.h5
    ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/C_diamond-twist-third-cart.h5
    fname Na2.structure.xml simple.txt)

foreach(fname cn.wfnoj.xml cn.wfj.xml cn.msd-wfnoj.xml cn.msd-1det-wfnoj.xml cn.msd-1det-wfj.xml cn.msd-wfj.xml C_diamond-twist-third-cart.wfj.xml)
  list(APPEND UTEST_INPUT_FILES ${qmcpack_SOURCE_DIR}/src/QMCWaveFunctions/tests/${fname})
endforeach()

foreach(fname C.BFD.xml Na.BFD.xml so_ecp_test.xml C.ccECP.xml N.ccECP.xml Be.BFD.xml)
  list(APPEND UTEST_INPUT_FILES ${qmcpack_SOURCE_DIR}/tests/pseudopotentials_for_tests/${fname})
endforeach()

foreach(CATEGORY coulomb force ham ewald2d)
  set(UTEST_EXE test_${SRC_DIR}_${CATEGORY})
  set(UTEST_NAME deterministic-unit_${UTEST_EXE})
  string(TOUPPER "${CATEGORY}_SRCS" SOURCE_FILE_VAR_NAME)
  add_executable(${UTEST_EXE} ${${SOURCE_FILE_VAR_NAME}})

  target_link_libraries(${UTEST_EXE} catch_main qmcham utilities_for_test minimal_test_pools)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} qmcwfs qmcwfs_omptarget qmcparticle qmcparticle_omptarget spline2 spline2_omptarget qmcutil platform_omptarget_LA)
  endif()

  add_unit_test(BASENAME ${UTEST_NAME} TEST_BINARY $<TARGET_FILE:${UTEST_EXE}> INPUT_FILES ${UTEST_INPUT_FILES})
endforeach()

# Run each fatal input independently; an abort must not terminate the regular Catch suite.
foreach(scenario missing_charge charge_above_tolerance charge_below_tolerance missing_atomic_number
                 atomic_number_mismatch atomic_number_alias_mismatch unknown_species unsupported_format factory)
  foreach(ranks 1 2)
    if(ranks GREATER 1 AND (NOT HAVE_MPI OR NOT scenario STREQUAL "factory"))
      continue()
    endif()
    set(abort_test deterministic-unit_test_ecp_builder_${scenario}-r${ranks}-t1)
    set(abort_dir "${CMAKE_CURRENT_BINARY_DIR}/${abort_test}")
    file(MAKE_DIRECTORY "${abort_dir}")
    configure_file("${qmcpack_SOURCE_DIR}/tests/pseudopotentials_for_tests/C.BFD.xml"
                   "${abort_dir}/C.BFD.xml" COPYONLY)
    set(abort_command ${QMC_GPU_TEST_LAUNCHER})
    if(HAVE_MPI)
      list(APPEND abort_command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND abort_command $<TARGET_FILE:test_hamiltonian_ham>)
    add_test(NAME ${abort_test}
             COMMAND ${CMAKE_COMMAND} "-DTEST_COMMAND=${abort_command}" "-DSCENARIO=${scenario}"
                     -P "${CMAKE_CURRENT_SOURCE_DIR}/check_ecp_builder_abort.cmake")
    set_tests_properties(${abort_test} PROPERTIES WORKING_DIRECTORY "${abort_dir}" TIMEOUT 45
                         PROCESSORS ${ranks} ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "unit")
    set_test_gpu_resources(${abort_test})
  endforeach()
endforeach()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
