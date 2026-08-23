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

list(APPEND SHORT_BCCH333_SCALARS "totenergy" "-28.88967350 0.1")
list(APPEND SHORT_BCCH333_SCALARS "kinetic" "38.55257604 0.5")
list(APPEND SHORT_BCCH333_SCALARS "potential" "-67.97693634 0.5")
list(APPEND SHORT_BCCH333_SCALARS "eeenergy" "-22.90804712 0.1")
list(APPEND SHORT_BCCH333_SCALARS "ionion" "-36.93144136 0.000001")
list(APPEND SHORT_BCCH333_SCALARS "samples" "1600 0.000000")

qmc_run_and_check(
  short-bccH_3x3x3_ae-vmc_bf
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
  qmc-short-bf
  qmc-short-bf.in.xml
  1
  4
  TRUE
  0
  SHORT_BCCH333_SCALARS # VMC
)

qmc_run_and_check(
  short-bccH_3x3x3_ae-vmcbatch_bf
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
  qmc-short-bf-batch
  qmc-short-bf-batch.in.xml
  1
  4
  TRUE
  0
  SHORT_BCCH333_SCALARS # VMC
)

list(APPEND DET_BCCH333_SCALARS "totenergy" "-29.43672031 0.5")
list(APPEND DET_BCCH333_SCALARS "kinetic" "33.78052518 2.5")
list(APPEND DET_BCCH333_SCALARS "potential" "-63.21724549 2.5")
list(APPEND DET_BCCH333_SCALARS "eeenergy" "-22.83921666 1.0")
list(APPEND DET_BCCH333_SCALARS "ionion" "-36.93144136 0.000001")
list(APPEND DET_BCCH333_SCALARS "samples" "20 0.000000")

qmc_run_and_check(
  deterministic-bccH_3x3x3_ae-vmc_bf
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
  det-qmc-short-bf
  det-qmc-bf.in.xml
  1
  1
  TRUE
  0
  DET_BCCH333_SCALARS # VMC
)

qmc_run_and_check(
  deterministic-bccH_3x3x3_ae-vmcbatch_bf
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
  det-qmc-short-bf-batch
  det-qmc-bf-batch.in.xml
  1
  1
  TRUE
  0
  DET_BCCH333_SCALARS # VMC
)
