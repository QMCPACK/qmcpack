#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//////////////////////////////////////////////////////////////////////////////////////
# Real-checkpoint smoke test; this checks execution and finite output, not accuracy.
set(_deepqmc_workdir "${CMAKE_CURRENT_BINARY_DIR}/deepqmc-he-vmc")
file(MAKE_DIRECTORY "${_deepqmc_workdir}")
configure_file("${qmcpack_SOURCE_DIR}/examples/molecules/He/he_deepqmc.xml" "${_deepqmc_workdir}/he_deepqmc.xml" COPYONLY)
maybe_symlink("${QMC_DEEPQMC_HE_CHECKPOINT}" "${_deepqmc_workdir}/chkpt.pt")
run_qmc_app_no_copy(deepqmc-he-vmc "${_deepqmc_workdir}" 1 1 _deepqmc_test_added he_deepqmc.xml)
if(_deepqmc_test_added)
  deepqmc_set_test_environment(deepqmc-he-vmc)
  set_property(TEST deepqmc-he-vmc APPEND PROPERTY LABELS "deepqmc;integration")
  set_tests_properties(deepqmc-he-vmc PROPERTIES TIMEOUT 120 FIXTURES_SETUP deepqmc_he_run)
  add_test(NAME deepqmc-he-vmc-check
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_scalar.py" "${_deepqmc_workdir}")
  set_tests_properties(deepqmc-he-vmc-check PROPERTIES
    LABELS "deepqmc;integration" TIMEOUT 120 FIXTURES_REQUIRED deepqmc_he_run)
endif()
