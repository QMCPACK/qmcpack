MESSAGE("Adding tests for afqmctools and QMCPACK workflow")

INCLUDE("${qmcpack_SOURCE_DIR}/CMake/test_labels.cmake")

FUNCTION(ADD_AFQMC_PYSCF_WORKFLOW_TEST test_name input_dir TEST_SCALARS TEST_OBS NPROC_INT)
  SET(TESTDIR "${CMAKE_CURRENT_BINARY_DIR}/converter_test_qmcpack_${test_name}-16-1")
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory "${input_dir}" "${TESTDIR}")

  # 1. Run scf.
  SET(TESTNAME "converter_test_pyscf_${test_name}")
  ADD_TEST(NAME ${TESTNAME} COMMAND ${TESTDIR}/scf.py gen_test_data WORKING_DIRECTORY ${TESTDIR})
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;afqmc")

  # 2.a Generate hamiltonian/wavefunction.
  SET(INT_SCRIPT)
  SET(TESTNAME "converter_test_pyscf_to_afqmc_${test_name}")
  ADD_TEST(NAME ${TESTNAME}
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC_INT} ${MPIEXEC_PREFLAGS}
           ${PROJECT_SOURCE_DIR}/utils/afqmctools/bin/pyscf_to_afqmc.py ${ARGN}
           WORKING_DIRECTORY ${TESTDIR})
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "converter_test_pyscf_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;afqmc" )

  # 2.b Test Hamiltonian consistency.
  SET(TESTNAME "converter_test_hamiltonian_consistency_${test_name}")
  ADD_TEST(NAME ${TESTNAME} COMMAND
           ${PROJECT_SOURCE_DIR}/utils/afqmctools/bin/test_afqmc_input.py --input afqmc.h5
           WORKING_DIRECTORY ${TESTDIR})
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "converter_test_pyscf_to_afqmc_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;afqmc")

  # 3. Generate xml input.
  SET(TESTNAME "converter_test_gen_input_${test_name}")
  ADD_TEST(NAME ${TESTNAME} COMMAND "${TESTDIR}/gen_input.py"
           WORKING_DIRECTORY ${TESTDIR})
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "converter_test_hamiltonian_consistency_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;afqmc")

  # 5. Run AFQMC.
  SET(TESTNAME "converter_test_qmcpack_${test_name}")
  QMC_RUN_AND_CHECK(${TESTNAME}
                    ${input_dir}
                    qmc afqmc.xml
                    16 1
                    TRUE
                    0 TEST_SCALARS)
  SET_TESTS_PROPERTIES(${TESTNAME}-16-1 PROPERTIES  DEPENDS "converter_test_pyscf_to_afqmc_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME}-16-1 APPEND PROPERTY LABELS "converter;afqmc")
  IF(TEST_OBS)
    SET(TESTNAME "converter_test_qmcpack_obs_${test_name}")
    ADD_TEST(NAME ${TESTNAME} COMMAND ${TESTDIR}/reference/test_obs.py
             WORKING_DIRECTORY ${TESTDIR})
    SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "converter_test_qmcpack_${test_name}-16-1")
    SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "converter;afqmc")
  ENDIF()
ENDFUNCTION()

SET(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples/afqmc)
# Molecular tests can be run through using either real or complex builds. No CUDA yet.
IF(NOT QMC_COMPLEX)
  SET(REAL_FLAG "-r")
ELSE()
  SET(REAL_FLAG "")
ENDIF()

IF(NOT ENABLE_CUDA)
  LIST(APPEND TEST_SCALARS_NEON "EnergyEstim__nume_real" "-128.71175210 0.00361907")
  ADD_AFQMC_PYSCF_WORKFLOW_TEST("01-neon_atom"
                                "${EXAMPLE_DIR}/01-neon_atom"
                                "${TEST_SCALARS_NEON}" FALSE 1
                                -v -t 1e-5 ${REAL_FLAG} -i reference/scf_ref.chk -o afqmc.h5)
  LIST(APPEND TEST_SCALARS_NEON "EnergyEstim__nume_real" "-128.70956768 0.00469893")
  ADD_AFQMC_PYSCF_WORKFLOW_TEST("02-neon_frozen_core"
                                "${EXAMPLE_DIR}/02-neon_frozen_core"
                                "${TEST_SCALARS_NEON}" FALSE 1
                                -v -t 1e-5 ${REAL_FLAG} -i reference/scf_ref.chk -o afqmc.h5 -c 8,22)
  LIST(APPEND TEST_SCALARS_C "EnergyEstim__nume_real" " -37.78271861 0.00287917")
  ADD_AFQMC_PYSCF_WORKFLOW_TEST("03-carbon_triplet_uhf"
                                "${EXAMPLE_DIR}/03-carbon_triplet_uhf"
                                "${TEST_SCALARS_C}" FALSE 1
                                -v -t 1e-5 ${REAL_FLAG} -i reference/scf_ref.chk -a -o afqmc.h5)
  IF(QMC_COMPLEX)
      # Avoid pyscf_to_afqmc to generate the hamiltonian for examples 4 and 5.
      # Due to test setup can't test real build.
      LIST(APPEND TEST_SCALARS_N2_NOMSD "EnergyEstim__nume_real" "-109.01813744 0.01562922")
      ADD_AFQMC_PYSCF_WORKFLOW_TEST("04-N2_nomsd"
                                    "${EXAMPLE_DIR}/04-N2_nomsd"
                                    "${TEST_SCALARS_N2_NOMSD}" FALSE 1 -h)
      LIST(APPEND TEST_SCALARS_N2_PHMSD "EnergyEstim__nume_real" "-109.04700629 0.01610529")
      ADD_AFQMC_PYSCF_WORKFLOW_TEST("05-N2_phmsd"
                                    "${EXAMPLE_DIR}/05-N2_phmsd"
                                    "${TEST_SCALARS_N2_PHMSD}" FALSE 1 -h)
      LIST(APPEND TEST_SCALARS_BP "EnergyEstim__nume_real" "-39.80583013 0.00134888")
      ADD_AFQMC_PYSCF_WORKFLOW_TEST("06-methane_converge_back_prop"
                                    "${EXAMPLE_DIR}/06-methane_converge_back_prop"
                                    "${TEST_SCALARS_BP}" TRUE 1
                                    -v -t 1e-5 -i reference/scf_ref.chk -o afqmc.h5)
      LIST(APPEND TEST_SCALARS_SC "EnergyEstim__nume_real" "-88.43097328 0.04224239")
      ADD_AFQMC_PYSCF_WORKFLOW_TEST("07-diamond_2x2x2_supercell"
                                    "${EXAMPLE_DIR}/07-diamond_2x2x2_supercell"
                                    "${TEST_SCALARS_SC}" FALSE 8
                                    -v -t 1e-5 -i reference/scf_ref.chk -o afqmc.h5)
      # Only kpoint-symmetric factorization has cuda build.
      LIST(APPEND TEST_SCALARS_SC "EnergyEstim__nume_real" "-88.44303023 0.04425413")
      ADD_AFQMC_PYSCF_WORKFLOW_TEST("08-diamond_2x2x2_kpoint_sym"
                                    "${EXAMPLE_DIR}/08-diamond_2x2x2_kpoint_sym"
                                    "${TEST_SCALARS_SC}" FALSE 8
                                    -v -t 1e-5 -k -a -i reference/scf_ref.chk -o afqmc.h5)
  ENDIF()
ENDIF()
