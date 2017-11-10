
# Add tests for Nexus examples

INCLUDE("${qmcpack_SOURCE_DIR}/CMake/python.cmake")

FUNCTION(ADD_EXAMPLE_TEST test_name test_dir nexus_file)
    #MESSAGE("Adding Nexus test ${test_name}")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}" "${CMAKE_CURRENT_BINARY_DIR}/${test_dir}")
    ADD_TEST(NAME nexus_example_${test_name} COMMAND python nexus_example_test.py "${CMAKE_CURRENT_BINARY_DIR}/${test_dir}" --nexus-file "${nexus_file}" --nexus-path "${CMAKE_CURRENT_SOURCE_DIR}/../library")
    SET_TESTS_PROPERTIES(nexus_example_${test_name} PROPERTIES TIMEOUT 240)
    SET_TESTS_PROPERTIES(nexus_example_${test_name} PROPERTIES LABELS "nexus")

ENDFUNCTION()

SET(NEXUS_PYTHON_REQS numpy)
CHECK_PYTHON_REQS(NEXUS_PYTHON_REQS "" ADD_TEST)

if (ADD_TEST)
  MESSAGE("Adding Nexus example tests")
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/nexus_example_test.py" "${CMAKE_CURRENT_BINARY_DIR}")
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/qmcpack/pseudopotentials" "${CMAKE_CURRENT_BINARY_DIR}/qmcpack/pseudopotentials")
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/quantum_espresso/pseudopotentials" "${CMAKE_CURRENT_BINARY_DIR}/quantum_espresso/pseudopotentials")
  ADD_EXAMPLE_TEST(c20 qmcpack/c20 c20.py)
  ADD_EXAMPLE_TEST(diamond qmcpack/diamond diamond.py)
  ADD_EXAMPLE_TEST(diamond_vacancy qmcpack/diamond diamond_vacancy.py)
  ADD_EXAMPLE_TEST(graphene qmcpack/graphene graphene.py)
  ADD_EXAMPLE_TEST(H2O qmcpack/H2O H2O.py)
  ADD_EXAMPLE_TEST(LiH qmcpack/LiH LiH.py)
  ADD_EXAMPLE_TEST(oxygen_dimer qmcpack/oxygen_dimer oxygen_dimer.py)
  ADD_EXAMPLE_TEST(relax_Ge_T_vs_kpoints    quantum_espresso/relax_Ge_T_vs_kpoints relax_vs_kpoints_example.py)
ELSE()
  MESSAGE("Skipping Nexus example tests because numpy is not present in python installation")
ENDIF()
