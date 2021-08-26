# Deterministic test for QMC optimizer

if(MIXED_PRECISION)
  message(VERBOSE "Skipping optimizer deterministic He tests in mixed precision (QMC_MIXED_PRECISION=1)")
elseif(QMC_CUDA)
  message(VERBOSE "Skipping optimizer deterministic He tests because CUDA is not supported (QMC_CUDA=1)")
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

  list(APPEND He_VMC_SCALARS "totenergy" "-3.25291030 0.000001")

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

  # Test batched version of optimizer

  list(APPEND He_batch_VMC_SCALARS "totenergy" "-2.2273099 0.000001")

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

endif()
