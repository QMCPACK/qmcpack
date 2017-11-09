
IF(NOT MIXED_PRECISION)
  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-grad-lap
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    grad_lap.xml
    1 1
    ../short-bccH_3x3x3_ae-deriv-1-1/check_grad_lap.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    deriv.xml
    1 1
    ../short-bccH_3x3x3_ae-deriv-1-1/check_deriv.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-gamma-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    gamma_deriv.xml
    1 1
    ../short-bccH_3x3x3_ae-deriv-1-1/check_deriv.py
  )
ENDIF(NOT MIXED_PRECISION)
