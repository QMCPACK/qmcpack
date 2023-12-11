if(NOT QMC_COMPLEX)

  # VMC short run with no Jastrows and UN-Truncated WF
  # Compare directly to sCI Variational Energy of: E=-8.04250444504137
  #Data generated using utils/make_ref_data.sh  LiH_ae_MSD_NOJ_SHORT qmc-ref/vmc_short_noj.s000.scalar.dat  qmc-ref/qmc-ref-vmc-long-noj.s000.scalar.dat
  list(APPEND LiH_ae_MSD_NOJ_SHORT_SCALARS "totenergy" "-8.042553312461 0.002407485304")
  list(APPEND LiH_ae_MSD_NOJ_SHORT_SCALARS "variance" "1.245679324458 0.051278600983")
  list(APPEND LiH_ae_MSD_NOJ_SHORT_SCALARS "kinetic" "8.004096666221 0.066070685701")
  list(APPEND LiH_ae_MSD_NOJ_SHORT_SCALARS "potential" "-16.046649978685 0.066254492826")
  list(APPEND LiH_ae_MSD_NOJ_SHORT_SCALARS "eeenergy" "3.379729221478 0.003550933783")
  list(APPEND LiH_ae_MSD_NOJ_SHORT_SCALARS "ionion" "0.995380119260 0.001000000000")
  qmc_run_and_check(
    short-LiH_ae-vmc_msdj_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    vmc_short_noj
    vmc_short_noj.in.xml
    1
    16
    TRUE
    0
    LiH_ae_MSD_NOJ_SHORT_SCALARS #VMC
  )

  #
  # VMC short run with Jastrows and Truncated WF
  # Reference VMC run in qmc-ref-vmc-short.s000.scalar.dat
  # VMC short run with Jastrow with 184 Determinants out  of 1847.
  #Data generated using utils/make_ref_data.sh  LiH_ae_MSD qmc-ref/vmc_short.s000.scalar.dat qmc-ref/qmc-ref-vmc-long.s000.scalar.dat
  list(APPEND LiH_ae_MSD_SHORT_SCALARS "totenergy" "-8.035845817205 0.000819929165")
  list(APPEND LiH_ae_MSD_SHORT_SCALARS "variance" "0.255494207230 0.003183663838")
  list(APPEND LiH_ae_MSD_SHORT_SCALARS "kinetic" "8.151755889083 0.036091574258")
  list(APPEND LiH_ae_MSD_SHORT_SCALARS "potential" "-16.187601706300 0.035811850526")
  list(APPEND LiH_ae_MSD_SHORT_SCALARS "eeenergy" "3.262235428602 0.002235254233")
  list(APPEND LiH_ae_MSD_SHORT_SCALARS "ionion" "0.995380119260 0.001000000000")
  qmc_run_and_check(
    short-LiH_ae-vmc_msdj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    vmc_short
    vmc_short.in.xml
    1
    16
    TRUE
    0
    LiH_ae_MSD_SHORT_SCALARS # VMC
  )

  #
  # VMC short run with Jastrows and optimized 184 determinants read from previously saved coefficients
  # Purpose: Reading and Replacing optimized Coeffs from previously optimized run.
  #Data generated using utils/make_ref_data.sh make_ref_data.sh LiH_ae_MSD_OPT_CI_SHORT  qmc-ref/vmc_short_opt_CI.s000.scalar.dat  qmc-ref/qmc-ref-vmc-long-opt-CI.s000.scalar.dat
  list(APPEND LiH_ae_MSD_OPT_CI_SHORT_SCALARS "totenergy" "-8.065033555875 0.000853294542")
  list(APPEND LiH_ae_MSD_OPT_CI_SHORT_SCALARS "variance" "0.053562868983 0.003478707706")
  list(APPEND LiH_ae_MSD_OPT_CI_SHORT_SCALARS "kinetic" "8.099258390791 0.045167286002")
  list(APPEND LiH_ae_MSD_OPT_CI_SHORT_SCALARS "potential" "-16.164291946667 0.044848012272")
  list(APPEND LiH_ae_MSD_OPT_CI_SHORT_SCALARS "eeenergy" "3.346222009202 0.002397327879")
  list(APPEND LiH_ae_MSD_OPT_CI_SHORT_SCALARS "ionion" "0.995380119260 0.001000000000")
  qmc_run_and_check(
    short-LiH_ae-vmc_msdj_optCI
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    vmc_short_opt_CI
    vmc_short_opt_CI.in.xml
    1
    16
    TRUE
    0
    LiH_ae_MSD_OPT_CI_SHORT_SCALARS # VMC
  )

  #
  # VMC Long run with Jastrows and Truncated WF
  # Reference VMC run in qmc-ref-vmc-long.s000.scalar.dat
  # VMC long run with Jastrow with 184 Determinants out  of 1847.
  #Data generated using utils/make_ref_data.sh  LiH_ae_MSD qmc-ref/vmc_long.s000.scalar.dat qmc-ref/qmc-ref-vmc-long.s000.scalar.dat
  list(APPEND LiH_ae_MSD_SCALARS "totenergy" "-8.035845817205 0.000266400021")
  list(APPEND LiH_ae_MSD_SCALARS "variance" "0.255494207230 0.000880077660")
  list(APPEND LiH_ae_MSD_SCALARS "kinetic" "8.151755889083 0.015127181063")
  list(APPEND LiH_ae_MSD_SCALARS "potential" "-16.187601706300 0.015063626138")
  list(APPEND LiH_ae_MSD_SCALARS "eeenergy" "3.262235428602 0.000802865851")
  list(APPEND LiH_ae_MSD_SCALARS "ionion" "0.995380119260 0.001000000000")
  qmc_run_and_check(
    long-LiH_ae-vmc_msdj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    vmc_long
    vmc_long.in.xml
    1
    16
    TRUE
    0
    LiH_ae_MSD_SCALARS # VMC
  )

  #
  # VMC Long run with no Jastrows and UN-Truncated WF
  # Compare directly to sCI Variational Energy of: E=-8.04250444504137
  #Data generated using utils/make_ref_data.sh  LiH_ae_MSD_NOJ qmc-ref/vmc_long_noj.s000.scalar.dat  qmc-ref/qmc-ref-vmc-long-noj.s000.scalar.dat
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "HartreeFock" "-8.04250444504137 0.0")
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "totenergy" "-8.042553312461 0.000488533041")
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "variance" "1.245679324458 0.048598851230")
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "kinetic" "8.004096666221 0.014528024773")
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "potential" "-16.046649978685 0.014554358540")
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "eeenergy" "3.379729221478 0.000841913994")
  list(APPEND LiH_ae_MSD_NOJ_SCALARS "ionion" "0.995380119260 0.001000000000")
  qmc_run_and_check(
    long-LiH_ae-vmc_msdj_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    vmc_long_noj
    vmc_long_noj.in.xml
    1
    16
    TRUE
    0
    LiH_ae_MSD_NOJ_SCALARS # VMC
  )
  #
  # VMC Long run with Jastrows and optimized 184 determinants read from previously saved coefficients
  # Purpose: Reading and Replacing optimized Coeffs from previously optimized run.
  #Data generated using utils/make_ref_data.sh make_ref_data.sh LiH_ae_MSD_OPT_CI  qmc-ref/vmc_long_opt_CI.s000.scalar.dat  qmc-ref/qmc-ref-vmc-long-opt-CI.s000.scalar.dat
  list(APPEND LiH_ae_MSD_OPT_CI_SCALARS "totenergy" "-8.065033555875 0.000289417316")
  list(APPEND LiH_ae_MSD_OPT_CI_SCALARS "variance" "0.053562868983 0.000909787243")
  list(APPEND LiH_ae_MSD_OPT_CI_SCALARS "kinetic" "8.099258390791 0.012893891746")
  list(APPEND LiH_ae_MSD_OPT_CI_SCALARS "potential" "-16.164291946667 0.012803640465")
  list(APPEND LiH_ae_MSD_OPT_CI_SCALARS "eeenergy" "3.346222009202 0.000692976585")
  list(APPEND LiH_ae_MSD_OPT_CI_SCALARS "ionion" "0.995380119260 0.001000000000")
  qmc_run_and_check(
    long-LiH_ae-vmc_msdj_optCI
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    vmc_long_opt_CI
    vmc_long_opt_CI.in.xml
    1
    16
    TRUE
    0
    LiH_ae_MSD_OPT_CI_SCALARS # VMC
  )

  # Minimal test of ensemble run
  # No check on results
  # Arrange for "imbalanced" ensemble and assume will crash or hang/timeout on major error
  # Should abort if number of MPI tasks is smaller than ensemble
  qmc_run_and_check(
    deterministic-ensemble_blocks_LiH_ae-vmc_msdj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    det_ensemble_blocks_vmc_msdj
    det_ensemble_blocks_vmc_msdj.txt
    4
    1
    TRUE)

  qmc_run_and_check(
    deterministic-ensemble_blocks_LiH_ae-vmc_msdj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    det_ensemble_blocks_vmc_msdj
    det_ensemble_blocks_vmc_msdj.txt
    3
    1
    FALSE)

  qmc_run_and_check(
    deterministic-ensemble_blocks_LiH_ae-vmc_msdj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD"
    det_ensemble_blocks_vmc_msdj
    det_ensemble_blocks_vmc_msdj.txt
    5
    1
    TRUE)

else()
  message(VERBOSE
          "Skipping LiH_ae_MSD tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
