if(NOT MIXED_PRECISION)
  simple_run_and_check(short-bccH_3x3x3_ae-deriv "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae" deriv.xml 1 1
                       check_deriv.py)

  simple_run_and_check(short-bccH_3x3x3_ae-grad_lap "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae" grad_lap.xml 1 1
                       check_grad_lap.py)

  simple_run_and_check(short-bccH_3x3x3_ae-gamma-deriv "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
                       gamma_deriv.xml 1 1 check_deriv.py)

  simple_run_and_check(short-bccH_3x3x3_ae-not_orth-deriv "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
                       not-orth_deriv.xml 1 1 check_deriv.py)
else()
  message(VERBOSE "Skipping bccH_3x3x3_ae derivative tests in mixed precision (QMC_MIXED_PRECISION=1)")
endif(NOT MIXED_PRECISION)
