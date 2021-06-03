
IF (NOT MIXED_PRECISION)
  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-grad_lap
    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    grad_lap.xml
    1 1
    check_grad_lap.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-deriv
    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    deriv.xml
    1 1
    check_deriv.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-gamma-deriv
    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    gamma_deriv.xml
    1 1
    check_deriv.py
  )
ELSE()
 MESSAGE("Skipping bccH_2x2x2_ae derivative tests in mixed precision (QMC_MIXED_PRECISION=1)")
ENDIF (NOT MIXED_PRECISION)
