if(NOT QMC_CUDA)
  if(NOT MIXED_PRECISION)

  if(NOT QMC_COMPLEX)
   set(OFEXT "_real")
  else()
   set(OFEXT "_comp")
  endif()

  simple_run_and_check(short-bccH_2x2x2_ae-grad_lap "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae" grad_lap.xml 1 1
                       check_grad_lap.py)

  simple_run_and_check(short-bccH_2x2x2_ae-deriv "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae" deriv.xml 1 1
                       check_deriv.py)

  simple_run_and_check(short-bccH_2x2x2_ae-gamma-deriv "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
                       gamma_deriv.xml 1 1 check_deriv.py)

  simple_run_and_check(
    deterministic-bccH_2x2x2_ae_bf-opt
    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    det_qmc-bf-opt.in.xml
    1
    1
    check_opt.py
    -s
    1
    -p
    det_qmc-bf-opt
    -r
    ./qmc-ref/ref_det_qmc-bf-opt${OFEXT}.s001.opt.xml
    )

  list(APPEND SHORT_BCCH_SCALARS "totenergy" "-6.60791746 0.0163198")
  list(APPEND SHORT_BCCH_SCALARS "kinetic" "11.38644428 0.0049810")
  list(APPEND SHORT_BCCH_SCALARS "potential" "-17.99436173 0.017268")
  list(APPEND SHORT_BCCH_SCALARS "eeenergy" "-6.20440454 0.0094955")
  list(APPEND SHORT_BCCH_SCALARS "ionion" "-10.94265693 0.000001")

  qmc_run_and_check(
    short-bccH_2x2x2_ae-vmc_bf
    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    qmc-short-bf
    qmc-short-bf.in.xml
    1
    16
    TRUE
    0
    SHORT_BCCH_SCALARS # VMC
    )


  IF (QMC_COMPLEX)
   list(APPEND DET_BCCH_SCALARS "totenergy" "-4.58593001 0.000005")
   list(APPEND DET_BCCH_SCALARS "kinetic" "12.70085823 0.000005")
   list(APPEND DET_BCCH_SCALARS "potential" "-17.28678824 0.000005")
   list(APPEND DET_BCCH_SCALARS "eeenergy" "-6.96009998 0.000005")
   list(APPEND DET_BCCH_SCALARS "ionion" "-10.94265693 0.000005")
   list(APPEND DET_BCCH_SCALARS "samples" "9 0.0")
  ELSE()
   list(APPEND DET_BCCH_SCALARS "totenergy" "-7.24345470 0.000005")
   list(APPEND DET_BCCH_SCALARS "kinetic" "11.08293632 0.000005")
   list(APPEND DET_BCCH_SCALARS "potential" "-18.32639101 0.000005")
   list(APPEND DET_BCCH_SCALARS "eeenergy" "-6.98494626 0.000005")
   list(APPEND DET_BCCH_SCALARS "ionion" "-10.94265693 0.000005")
   list(APPEND DET_BCCH_SCALARS "samples" "9 0.0")
  ENDIF()

  qmc_run_and_check(
    deterministic-bccH_2x2x2_ae-vmc_bf
    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    det_qmc-bf-opt
    det_qmc-bf-opt.in.xml
    1
    1
    TRUE
    2
    DET_BCCH_SCALARS # VMC
    )

  else()
    message_verbose("Skipping bccH_2x2x2_ae derivative and optimization tests in mixed precision (QMC_MIXED_PRECISION=1)")
  endif(NOT MIXED_PRECISION)
else()
  message_verbose("Skipping bccH_2x2x2_ae derivative and backflow tests because they are not support by CUDA build (QMC_CUDA=1)")
endif()

