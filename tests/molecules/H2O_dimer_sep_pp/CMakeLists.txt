
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

# 2 water molecules far away test for Tmoves.
# short tests
  LIST(APPEND H2OS2_DMCLA_SCALARS "totenergy" "-34.550288  0.0015") # total energy
  LIST(APPEND H2OS2_DMCLA_SCALARS "nonlocalecp" "3.7922 0.012") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_la
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-LA.xml
                    1 16
                    ${MP_SUCCESS}
                    2 H2OS2_DMCLA_SCALARS
                    )

  LIST(APPEND H2OS2_DMCTM0_SCALARS "totenergy" "-34.556408  0.0015") # total energy
  LIST(APPEND H2OS2_DMCTM0_SCALARS "nonlocalecp" "3.3476 0.012") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_tm0
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-TM0.xml
                    1 16
                    ${MP_SUCCESS}
                    2 H2OS2_DMCTM0_SCALARS
                    )

  LIST(APPEND H2OS2_DMCTM1_SCALARS "totenergy" "-34.561071  0.0015") # total energy
  LIST(APPEND H2OS2_DMCTM1_SCALARS "nonlocalecp" "3.1809 0.012") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_tm1
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-TM1.xml
                    1 16
                    ${MP_SUCCESS}
                    2 H2OS2_DMCTM1_SCALARS
                    )

# long tests
  LIST(APPEND LONG_H2OS2_DMCLA_SCALARS "totenergy" "-34.550288  0.0005") # total energy
  LIST(APPEND LONG_H2OS2_DMCLA_SCALARS "nonlocalecp" "3.7922 0.0045") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_la
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-LA.xml
                    1 16
                    ${MP_SUCCESS}
                    2 LONG_H2OS2_DMCLA_SCALARS
                    )

  LIST(APPEND LONG_H2OS2_DMCTM0_SCALARS "totenergy" "-34.556408  0.0005") # total energy
  LIST(APPEND LONG_H2OS2_DMCTM0_SCALARS "nonlocalecp" "3.3476 0.0045") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_tm0
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-TM0.xml
                    1 16
                    ${MP_SUCCESS}
                    2 LONG_H2OS2_DMCTM0_SCALARS
                    )

  LIST(APPEND LONG_H2OS2_DMCTM1_SCALARS "totenergy" "-34.561071  0.0005") # total energy
  LIST(APPEND LONG_H2OS2_DMCTM1_SCALARS "nonlocalecp" "3.1809 0.0045") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_tm1
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-TM1.xml
                    1 16
                    ${MP_SUCCESS}
                    2 LONG_H2OS2_DMCTM1_SCALARS
                    )

 ELSE()
    MESSAGE("Skipping H2O_dimer_sep_pp tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping H2O_dimer_sep_pp tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

