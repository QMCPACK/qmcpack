
 IF (NOT QMC_CUDA)

# Reference VMC run in qmc-ref "-10.495941  0.000073"
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.495941 0.0065")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.51198 0.078")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.00792 0.079")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.68548 0.018")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.2011 0.090")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.65437 0.024")
  LIST(APPEND DIAMOND_SCALARS "samples" "16000 0.0")



  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )



#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.495941 0.0017")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "160000 0.0")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.532497 0.000470"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.532497 0.006127")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

# Deterministic test
  IF (QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.26918214 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "10.15418453 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-20.42336667 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.32915939 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77567016 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-6.76235562 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "0.44381850 0.000005")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
  ELSE()
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.26902520 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "10.15433651 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-20.42336170 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.32916126 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566743 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-6.76235283 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "0.44381981 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
  ENDIF()


  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    det_qmc_short
                    det_qmc_short.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_SCALARS # VMC
                    )


  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
