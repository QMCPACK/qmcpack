MESSAGE("Adding tests for afqmctools and QMCPACK workflow")

INCLUDE("${qmcpack_SOURCE_DIR}/CMake/test_labels.cmake")

FUNCTION(ADD_AFQMC_PYSCF_WORKFLOW_TEST test_name input_dir TEST_SCALARS)
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory "${input_dir}" "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  # 1. Run scf.
  SET(TESTNAME "pyscf_${test_name}")
  message(${CMAKE_CURRENT_BINARY_DIR})
  ADD_TEST(NAME ${TESTNAME} COMMAND python ${CMAKE_CURRENT_BINARY_DIR}/${test_name}/scf.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "workflow;afqmc")
  # 2. Generate hamiltonian/wavefunction.
  SET(INT_SCRIPT )
  SET(TESTNAME "pyscf_to_afqmc_${test_name}")
  ADD_TEST(NAME ${TESTNAME} COMMAND python ${PROJECT_SOURCE_DIR}/utils/afqmctools/bin/pyscf_to_afqmc.py ${ARGN}
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "pyscf_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "workflow;afqmc" )
  # 3. Generate xml input.
  SET(TESTNAME "gen_input_${test_name}")
  ADD_TEST(NAME ${TESTNAME} COMMAND python "${CMAKE_CURRENT_BINARY_DIR}/${test_name}/gen_input.py"
           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "pyscf_to_afqmc_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "workflow;afqmc")
  SET(TESTNAME "copy_files_${test_name}")
  # 4 Need to copy files from ${test_dir} to directory created by QMC_RUN_AND_CHECK
  ADD_TEST(NAME ${TESTNAME} COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/copy_directory.sh
           ${CMAKE_CURRENT_BINARY_DIR}/${test_name} ${CMAKE_CURRENT_BINARY_DIR}/qmcpack_${test_name}-1-1)
  SET_TESTS_PROPERTIES(${TESTNAME} PROPERTIES  DEPENDS "gen_input_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME} APPEND PROPERTY LABELS "workflow;afqmc")
  # 5. Run AFQMC.
  SET(TESTNAME "qmcpack_${test_name}")
  #LIST(APPEND TEST_SCALARS_NEON_${test_name} "EnergyEstim__nume_real" "-128.702807 0.007712")
  QMC_RUN_AND_CHECK(${TESTNAME}
                    ${CMAKE_CURRENT_BINARY_DIR}/${test_name}
                    qmc afqmc.xml
                    1 1
                    TRUE
                    0 TEST_SCALARS)
  SET_TESTS_PROPERTIES(${TESTNAME}-1-1 PROPERTIES  DEPENDS "copy_files_${test_name}")
  SET_PROPERTY(TEST ${TESTNAME}-1-1 APPEND PROPERTY LABELS "workflow;afqmc")
ENDFUNCTION()

if(NOT ENABLE_CUDA)
  SET(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples/afqmc)
  #message("EXAMPLE: ${EXAMPLE_DIR} ${qmcpack_BINARY_DIR}")
  IF(QMC_COMPLEX)
    LIST(APPEND TEST_SCALARS_NEON "EnergyEstim__nume_real" "-128.702807 0.007712")
    ADD_AFQMC_PYSCF_WORKFLOW_TEST("01-neon_atom" "${EXAMPLE_DIR}/01-neon_atom" "${TEST_SCALARS_NEON}" -v -t 1e-5 -i scf.chk -o afqmc.h5)
    ADD_AFQMC_PYSCF_WORKFLOW_TEST("02-neon_frozen_core" "${EXAMPLE_DIR}/02-neon_frozen_core" "${TEST_SCALARS_NEON}" -v -t 1e-5 -i scf.chk -o afqmc.h5 -c 8,22)
    LIST(APPEND TEST_SCALARS_C "EnergyEstim__nume_real" "-37.78416183 0.00623919")
    ADD_AFQMC_PYSCF_WORKFLOW_TEST("03-carbon_triplet_uhf" "${EXAMPLE_DIR}/03-carbon_triplet_uhf" "${TEST_SCALARS_C}" -v -t 1e-5 -i scf.chk -a -o afqmc.h5)
    #ADD_AFQMC_PYSCF_WORKFLOW_TEST("04-N2_non_ortho_msd " "${EXAMPLE_DIR}/04-N2_non_orth_msd" TEST_SCALARS -v -t 1e-5 -i scf.chk -o afqmc.h5 -c 8,22)
    #ADD_AFQMC_PYSCF_WORKFLOW_TEST("05-N2_ph_msd" "${EXAMPLE_DIR}/05-N2_ph_msd" TEST_SCALARS -v -t 1e-5 -i scf.chk -o afqmc.h5 -c 8,22)
  ENDIF()
endif()
