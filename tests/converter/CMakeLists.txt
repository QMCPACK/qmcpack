include("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")

# Add tests for conversion of GAMESS to QMCPACK input (via convert4qmc)

function(ADD_CONVERTER_TEST test_name)
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}"
                          "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  set(EXE_NAME $<TARGET_FILE:convert4qmc>)
  if(HAVE_MPI)
    set(EXE_NAME ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                 $<TARGET_FILE:convert4qmc>)
    string(REPLACE ";" " " EXE_NAME "${EXE_NAME}")
  endif()
  add_test(NAME converter_${test_name}
           COMMAND "${CMAKE_CURRENT_BINARY_DIR}/converter_test.py" "${CMAKE_CURRENT_BINARY_DIR}/${test_name}" --exe
                   "${EXE_NAME}" --h5diff "${HDF5_DIFF_EXECUTABLE}")
  set_tests_properties(converter_${test_name} PROPERTIES 
                       FAIL_REGULAR_EXPRESSION "  FAIL"
                       PASS_REGULAR_EXPRESSION "  pass"
                       TIMEOUT 120 LABELS "converter;deterministic")
endfunction()

function(ADD_PWCONVERTER_TEST test_name)
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}"
                          "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  set(EXE_NAME $<TARGET_FILE:convertpw4qmc>)
  add_test(NAME converter_${test_name}
           COMMAND "${CMAKE_CURRENT_BINARY_DIR}/pwconverter_test.py" "${CMAKE_CURRENT_BINARY_DIR}/${test_name}" --exe
                   "${EXE_NAME}" --h5diff "${HDF5_DIFF_EXECUTABLE}")
  set_tests_properties(converter_${test_name} PROPERTIES 
                       FAIL_REGULAR_EXPRESSION "  FAIL"
                       PASS_REGULAR_EXPRESSION "  pass" 
                       TIMEOUT 120 LABELS "converter;deterministic")
endfunction()

execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/converter_test.py"
                        "${CMAKE_CURRENT_BINARY_DIR}")
add_converter_test(test_He_sto3g)
add_converter_test(test_Be_sto3g)
add_converter_test(test_Be_ccd)
add_converter_test(test_O_ext)
add_converter_test(test_C_sto3g)
add_converter_test(test_HCNp)
add_converter_test(test_aldet1)
add_converter_test(test_aldet5)
add_converter_test(test_aldet6)
if(HDF5_DIFF_EXECUTABLE)
  add_converter_test(test_HDF5_FeCO6)
  add_converter_test(test_HDF5_Be_ccd)
  add_converter_test(test_Bi_dirac)
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/pwconverter_test.py"
                          "${CMAKE_CURRENT_BINARY_DIR}")
  add_pwconverter_test(test_NaCl_qbox)
  add_pwconverter_test(test_O2_pwscf)
  set(H5_REF ${qmcpack_SOURCE_DIR}/tests/solids/monoO_noncollinear_1x1x1_pp/o2_45deg_spins.pwscf.h5)
  maybe_symlink(${H5_REF} ${CMAKE_CURRENT_BINARY_DIR}/test_O2_pwscf/gold.h5)
else()
  message("Skipping converter tests with HDF output because executable h5diff was not found")
endif()
add_converter_test(test_LiH_pyscf)
add_converter_test(test_LiH_pyscf_UHF)
add_converter_test(test_LiH_qp)
add_converter_test(test_C4_MSD_excited)
add_converter_test(test_diamond2_rmg)
add_converter_test(test_O_rmg)
