if(NOT MIXED_PRECISION)
  simple_run_and_check(deterministic-heg2d-4_ae-deriv "${qmcpack_SOURCE_DIR}/tests/heg/heg2d_4" deriv.xml 1 1
                       check_deriv.py)

  simple_run_and_check(deterministic-heg2d-4_ae-grad_lap "${qmcpack_SOURCE_DIR}/tests/heg/heg2d_4" grad_lap.xml 1 1
                       check_grad_lap.py)

else()
  message(VERBOSE "Skipping eg2d derivative tests in mixed precision (QMC_MIXED_PRECISION=1)")
endif(NOT MIXED_PRECISION)
