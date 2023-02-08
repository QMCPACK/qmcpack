# Deterministic test for QMC optimizer

if(MIXED_PRECISION)
  message(VERBOSE "Skipping optimizer deterministic He tests in mixed precision (QMC_MIXED_PRECISION=1)")
else()
  simple_run_and_check(
    deterministic-He_ae-opt
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt.xml
    1
    1
    check_opt.py
    -s
    1
    -p
    det_He_opt
    -r
    reference.det_He_opt.s001.opt.xml)

  simple_run_and_check(
    deterministic-He_ae-opt
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_mwalkers.xml
    1
    16
    check_opt.py
    -s
    1
    -p
    det_He_opt_mwalkers
    -r
    reference.det_He_opt_mwalkers_1-16.s001.opt.xml)

  simple_run_and_check(
    deterministic-He_ae-opt
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_mwalkers.xml
    2
    3
    check_opt.py
    -s
    1
    -p
    det_He_opt_mwalkers
    -r
    reference.det_He_opt_mwalkers_2-3.s001.opt.xml)

  list(APPEND He_VMC_SCALARS "totenergy" "-3.25291030 0.000001")

  list(APPEND He_VMC_MWALKERS_1-16_SCALAR "totenergy" "-2.84983901 0.000001")

  list(APPEND He_VMC_MWALKERS_2-3_SCALAR "totenergy" "-2.81743526 0.000001")

  qmc_run_and_check(
    deterministic-He_ae-opt_vmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt
    det_He_opt.xml
    1
    1
    TRUE
    2
    He_VMC_SCALARS # VMC
  )

  qmc_run_and_check(
    deterministic-He_ae-opt_vmc-mwalkers
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_mwalkers
    det_He_opt_mwalkers.xml
    1
    16
    TRUE
    2
    He_VMC_MWALKERS_1-16_SCALAR # VMC
  )

  qmc_run_and_check(
    deterministic-He_ae-opt_vmc-mwalkers
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_mwalkers
    det_He_opt_mwalkers.xml
    2
    3
    TRUE
    2
    He_VMC_MWALKERS_2-3_SCALAR # VMC
  )

  # Test batched version of optimizer

  list(APPEND He_batch_VMC_SCALARS "totenergy" "-3.25291030 0.000001")

  list(APPEND He_batch_VMC_MWALKERS_1-16_SCALARS "totenergy" "-2.82592127 0.000001")

  list(APPEND He_batch_VMC_MWALKERS_2-3_SCALARS "totenergy" "-2.92537421 0.000001")

  simple_run_and_check(
    deterministic-He_ae-opt-batch
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_batch.xml
    1
    1
    check_opt.py
    -s
    1
    -p
    det_He_opt_batch
    -r
    reference.det_He_opt.s001.opt.xml)

  simple_run_and_check(
    deterministic-He_ae-opt-batch-mwalkers
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_batch_mwalkers.xml
    1
    16
    check_opt.py
    -s
    1
    -p
    det_He_opt_batch_mwalkers
    -r
    reference.det_He_opt_batch_mwalkers_1-16.s001.opt.xml)

  simple_run_and_check(
    deterministic-He_ae-opt-batch-mwalkers
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_batch_mwalkers.xml
    2
    3
    check_opt.py
    -s
    1
    -p
    det_He_opt_batch_mwalkers
    -r
    reference.det_He_opt_batch_mwalkers_2-3.s001.opt.xml)

  qmc_run_and_check(
    deterministic-He_ae-opt_vmc-batch
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_batch
    det_He_opt_batch.xml
    1
    1
    TRUE
    2
    He_batch_VMC_SCALARS # VMC
  )

  qmc_run_and_check(
    deterministic-He_ae-opt_vmc-batch-mwalkers
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_batch_mwalkers
    det_He_opt_batch_mwalkers.xml
    1
    16
    TRUE
    2
    He_batch_VMC_MWALKERS_1-16_SCALARS # VMC
  )

  qmc_run_and_check(
    deterministic-He_ae-opt_vmc-batch-mwalkers
    "${qmcpack_SOURCE_DIR}/tests/molecules/He_ae"
    det_He_opt_batch_mwalkers
    det_He_opt_batch_mwalkers.xml
    2
    3
    TRUE
    2
    He_batch_VMC_MWALKERS_2-3_SCALARS # VMC
  )

endif()
