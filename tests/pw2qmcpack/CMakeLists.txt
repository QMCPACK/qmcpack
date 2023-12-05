#############################################################
# Add tests to ctest
#############################################################

message("Adding tests for pw.x, pw2qmcpack.x and QMCPACK workflow")

if(NOT TEST_MAX_PROCS)
  set(TEST_MAX_PROCS 16)
  message("TEST_MAX_PROCS was unset. Set to 16")
endif()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
  set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()
#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
else()
  include("${qmcpack_SOURCE_DIR}/CMake/run_qe.cmake")
  set(LAST_TEST_NAME "NONE")
  set(QE_TEST_NAME "NONE")

  run_qe_test(
    qe-LiH-unpolarized
    ${qmcpack_SOURCE_DIR}/examples/solids/dft-inputs
    16
    16
    1
    16
    1
    1
    LiH
    QE_TEST_NAME)

  #MESSAGE("${QE_TEST_NAME}")
  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME)

  # These values are slightly different from those in tests/systems
  # because the h5 is generated in different ways.
  list(APPEND LIH_GAMMA_SCALARS "kinetic" "7.46326 0.05")
  list(APPEND LIH_GAMMA_SCALARS "localecp" "-11.78311 0.05")
  list(APPEND LIH_GAMMA_SCALARS "samples" "1600000 0.0")

  set(THIS_TEST_NAME "${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj")
  qmc_run_and_check(
    ${THIS_TEST_NAME}
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-gamma
    hf_vmc_LiH-gamma.xml
    1
    16
    TRUE
    0
    LIH_GAMMA_SCALARS # VMC
  )

  if(${THIS_TEST_NAME})
    set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  endif()

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME)

  # These values are slightly different from those in tests/systems
  # because the h5 is generated in different ways.
  list(APPEND LIH_X_SCALARS "kinetic" "7.33899 0.05")
  list(APPEND LIH_X_SCALARS "localecp" "-11.37256 0.05")
  list(APPEND LIH_X_SCALARS "samples" "1600000 0.0")

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-x
    hf_vmc_LiH-x.xml
    1
    16
    TRUE
    0
    LIH_X_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  run_qe_test(
    qe-LiH-unpolarized
    ${qmcpack_SOURCE_DIR}/examples/solids/dft-inputs
    16
    16
    2
    16
    2
    2
    LiH
    QE_TEST_NAME)

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-gamma
    hf_vmc_LiH-gamma.xml
    1
    16
    TRUE
    0
    LIH_GAMMA_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-x
    hf_vmc_LiH-x.xml
    1
    16
    TRUE
    0
    LIH_X_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  # test case for no wf_collect
  run_qe_test(
    qe-LiH-unpolarized-no-collect
    ${qmcpack_SOURCE_DIR}/examples/solids/dft-inputs-no-collect
    16
    16
    16
    16
    1
    1
    LiH
    QE_TEST_NAME)

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-gamma
    hf_vmc_LiH-gamma.xml
    1
    16
    TRUE
    0
    LIH_GAMMA_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-x
    hf_vmc_LiH-x.xml
    1
    16
    TRUE
    0
    LIH_X_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  run_qe_test(
    qe-LiH-unpolarized-no-collect
    ${qmcpack_SOURCE_DIR}/examples/solids/dft-inputs-no-collect
    16
    12
    12
    16
    2
    2
    LiH
    QE_TEST_NAME)

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-gamma
    hf_vmc_LiH-gamma.xml
    1
    16
    TRUE
    0
    LIH_GAMMA_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
    hf_vmc_LiH-x
    hf_vmc_LiH-x.xml
    1
    16
    TRUE
    0
    LIH_X_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  # test case for spin polarized case
  run_qe_test(
    qe-LiH-polarized
    ${qmcpack_SOURCE_DIR}/examples/solids/dft-inputs-polarized
    16
    12
    6
    16
    4
    3
    LiH
    QE_TEST_NAME)

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/fake_nspin2"
    hf_vmc_LiH-gamma
    hf_vmc_LiH-gamma.xml
    1
    16
    TRUE
    0
    LIH_GAMMA_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/fake_nspin2"
    hf_vmc_LiH-x
    hf_vmc_LiH-x.xml
    1
    16
    TRUE
    0
    LIH_X_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  # test case for no wf_collect
  run_qe_test(
    qe-LiH-polarized-no-collect
    ${qmcpack_SOURCE_DIR}/examples/solids/dft-inputs-polarized-no-collect
    16
    12
    12
    16
    3
    3
    LiH
    QE_TEST_NAME)

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/fake_nspin2"
    hf_vmc_LiH-gamma
    hf_vmc_LiH-gamma.xml
    1
    16
    TRUE
    0
    LIH_GAMMA_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  softlink_h5(${QE_TEST_NAME} LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/fake_nspin2"
    hf_vmc_LiH-x
    hf_vmc_LiH-x.xml
    1
    16
    TRUE
    0
    LIH_X_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  # test case for O atom
  # from a long run
  # totenergy
  # LocalEnergy           =        -15.870035 +/-         0.000016
  # Kinetic               =          11.85233 +/-          0.00014
  # ElecElec              =         10.368642 +/-         0.000050
  # LocalECP              =         -40.15639 +/-          0.00020
  # NonLocalECP           =          2.065389 +/-         0.000083

  list(APPEND O_ATOM_SCALARS "totenergy" "-15.870035 0.01")
  list(APPEND O_ATOM_SCALARS "kinetic" "11.85233 0.1")
  list(APPEND O_ATOM_SCALARS "localecp" "-40.15639 0.1")
  list(APPEND O_ATOM_SCALARS "nonlocalecp" "2.065389 0.04")
  list(APPEND O_ATOM_SCALARS "eeenergy" "10.368642 0.02")
  list(APPEND O_ATOM_SCALARS "samples" "64000 0.0")

  run_qe_test(
    qe-atomO-polarized-no-collect
    ${qmcpack_SOURCE_DIR}/tests/solids/atomO_pp/dft-inputs
    16
    0
    16
    1
    0
    1
    atomO
    QE_TEST_NAME)

  softlink_h5(${QE_TEST_NAME} atomO_pp-vmc_sdj-1-16 atomO atomO.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-atomO_pp-vmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/atomO_pp"
    qmc_short
    qmc_short.in.xml
    1
    16
    TRUE
    0
    O_ATOM_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-atomO_pp-vmc_sdj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

  # pw2casino Kinetic 11.64963538
  # measured   Kinetic 11.64798 +/- 0.00017
  list(APPEND O_ATOM_NOJ_SCALARS "kinetic" "11.64798 0.1")

  softlink_h5(${QE_TEST_NAME} atomO_pp-vmc_noj-1-16 atomO atomO.pwscf.h5 LAST_TEST_NAME)

  qmc_run_and_check(
    ${QE_TEST_NAME}-atomO_pp-vmc_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/atomO_pp"
    qmc_short_noj
    qmc_short_noj.in.xml
    1
    16
    TRUE
    0
    O_ATOM_NOJ_SCALARS # VMC
  )

  set_tests_properties(${QE_TEST_NAME}-atomO_pp-vmc_noj-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})

endif()
