if(NOT QMC_COMPLEX)
  #
  # LiH molecule gaussian tests, all electron using Gamess
  # Energy from Gamess: E=-7.9873236457
  #
  #   "kinetic" "7.991319   0.000066") # kinetic energy
  #   "totenergy" "-7.9873131  0.0000088 ") # total energy
  #   "eeenergy" "3.488823  0.000010") # e-e energy
  #   "samples" "79134720000 0.0") # samples

  list(APPEND LIH_SHORT "kinetic" "7.991319 0.012") # kinetic energy
  list(APPEND LIH_SHORT "totenergy" "-7.9873131 0.0015") # total energy
  list(APPEND LIH_SHORT "eeenergy" "3.488823 0.00086") # e-e energy
  list(APPEND LIH_SHORT "samples" "9600000 0.0") # samples
  qmc_run_and_check(
    short-LiH_ae_gms-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_gms"
    vmc_short_noj
    vmc_short_noj.in.xml
    1
    16
    TRUE
    0
    LIH_SHORT # VMC
  )

  list(APPEND LIH_LONG "kinetic" "7.991319 0.0045") # kinetic energy
  list(APPEND LIH_LONG "totenergy" "-7.9873131 0.00052") # total energy
  list(APPEND LIH_LONG "eeenergy" "3.488823 0.00027") # e-e energy
  list(APPEND LIH_LONG "samples" "96000000 0.0") # samples
  qmc_run_and_check(
    long-LiH_ae_qp-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_gms"
    vmc_long_noj
    vmc_long_noj.in.xml
    1
    16
    TRUE
    0
    LIH_LONG # VMC
  )

  # estimator tests
  include("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
  check_python_reqs("numpy;h5py" diamond-estimator add_tests)
  if(add_tests)
    simple_run_and_check(
      short-LiH_ae-vmc_hf_noj_estimator_energydensity_voronoi
      "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_gms"
      vmc_noj_edens_vor_short.in.xml
      4
      4
      check_stats.py
      -s
      0
      -q
      'energydensity,EDvoronoi'
      -e
      2
      -c
      2
      -p
      vmc_noj_edens_vor_short
      -r
      qmc-ref/vmc_noj_edens_vor_short.s000.stat_ref_energydensity.dat)

    simple_run_and_check(
      long-LiH_ae-vmc_hf_noj_estimator_energydensity_voronoi
      "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_gms"
      vmc_noj_edens_vor_long.in.xml
      4
      4
      check_stats.py
      -s
      0
      -q
      'energydensity,EDvoronoi'
      -e
      2
      -c
      2
      -p
      vmc_noj_edens_vor_long
      -r
      qmc-ref/vmc_noj_edens_vor_long.s000.stat_ref_energydensity.dat)
  endif()

else()
  message(VERBOSE
          "Skipping LiH_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
