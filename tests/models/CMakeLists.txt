  LIST(APPEND SHO_SCALARS "totenergy" "0.45  0")
  LIST(APPEND SHO_SCALARS "kinetic"   "0.225 0")
  LIST(APPEND SHO_SCALARS "potential" "0.225 0")
  LIST(APPEND SHO_SCALARS "latdev"    "1.0   0")

  QMC_RUN_AND_CHECK(short-sho-vmc
    "${CMAKE_SOURCE_DIR}/tests/models/sho"
    m5
    vmc_sho.xml
    1 16
    TRUE
    0 BCC_H_SCALARS # VMC
  )
  
  SIMPLE_RUN_AND_CHECK(short-sho-grad_lap
    "${CMAKE_SOURCE_DIR}/tests/models/sho"
    grad_lap.xml
    1 1
    ../../solids/bccH_3x3x3_ae/short-bccH_3x3x3_ae-deriv-1-1/check_grad_lap.py
  )
