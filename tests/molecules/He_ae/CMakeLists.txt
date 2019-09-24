# Deterministic test for QMC optimizer 

IF (NOT MIXED_PRECISION)
  SIMPLE_RUN_AND_CHECK(
    deterministic-He_ae-opt
    "${CMAKE_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt.xml
    1 1
    check_opt.py -s 1 -p det_He_opt -r qmc-ref/det_He_opt.s001.opt.xml
    )
ELSE()
  MESSAGE_VERBOSE("Skipping optimizer deterministic tests in mixed precision (QMC_MIXED_PRECISION=1)")
ENDIF()

