
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

# 2 water molecules far away test for Tmoves.
# short tests
  LIST(APPEND H2OS2_DMCLA_SCALARS "totenergy" "-34.549539898120  0.001538707986") # total energy
  LIST(APPEND H2OS2_DMCLA_SCALARS "nonlocalecp" "3.796204397761 0.016487176409") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_la
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-LA.xml
                    1 16
                    TRUE
                    2 H2OS2_DMCLA_SCALARS
                    )

  LIST(APPEND H2OS2_DMCTM0_SCALARS "totenergy" "-34.556070141054  0.001522544760") # total energy
  LIST(APPEND H2OS2_DMCTM0_SCALARS "nonlocalecp" "3.341925103599 0.013031467636") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_tm0
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-TM0.xml
                    1 16
                    TRUE
                    2 H2OS2_DMCTM0_SCALARS
                    )

  LIST(APPEND H2OS2_DMCTM1_SCALARS "totenergy" "-34.560423250334  0.001457944721") # total energy
  LIST(APPEND H2OS2_DMCTM1_SCALARS "nonlocalecp" "3.175475067660 0.013122004395") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_tm1
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-TM1.xml
                    1 16
                    TRUE
                    2 H2OS2_DMCTM1_SCALARS
                    )

  LIST(APPEND H2OS2_DMCTM3_SCALARS "totenergy" "-34.560683215187  0.001479939483") # total energy
  LIST(APPEND H2OS2_DMCTM3_SCALARS "nonlocalecp" "3.183908083459 0.013476329953") # non local ecp

  QMC_RUN_AND_CHECK(short-H2O_dimer_sep_pp-j3_dmc_tm3
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    short-H2O-TM3.xml
                    1 16
                    TRUE
                    2 H2OS2_DMCTM3_SCALARS
                    )

# long tests
  LIST(APPEND LONG_H2OS2_DMCLA_SCALARS "totenergy" "-34.549539898120  0.000489581983") # total energy
  LIST(APPEND LONG_H2OS2_DMCLA_SCALARS "nonlocalecp" "3.796204397761 0.005245845604") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_la
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-LA.xml
                    1 16
                    TRUE
                    2 LONG_H2OS2_DMCLA_SCALARS
                    )

  LIST(APPEND LONG_H2OS2_DMCTM0_SCALARS "totenergy" "-34.556070141054  0.000484439211") # total energy
  LIST(APPEND LONG_H2OS2_DMCTM0_SCALARS "nonlocalecp" "3.341925103599 0.004146317448") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_tm0
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-TM0.xml
                    1 16
                    TRUE
                    2 LONG_H2OS2_DMCTM0_SCALARS
                    )

  LIST(APPEND LONG_H2OS2_DMCTM1_SCALARS "totenergy" "-34.560423250334  0.000463884944") # total energy
  LIST(APPEND LONG_H2OS2_DMCTM1_SCALARS "nonlocalecp" "3.175475067660 0.004175124192") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_tm1
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-TM1.xml
                    1 16
                    TRUE
                    2 LONG_H2OS2_DMCTM1_SCALARS
                    )

  LIST(APPEND LONG_H2OS2_DMCTM3_SCALARS "totenergy" "-34.560683215187  0.000470883178") # total energy
  LIST(APPEND LONG_H2OS2_DMCTM3_SCALARS "nonlocalecp" "3.183908083459 0.004287862548") # non local ecp

  QMC_RUN_AND_CHECK(long-H2O_dimer_sep_pp-j3_dmc_tm3
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
                    H2O
                    long-H2O-TM3.xml
                    1 16
                    TRUE
                    2 LONG_H2OS2_DMCTM3_SCALARS
                    )

 ELSE()
    MESSAGE("Skipping H2O_dimer_sep_pp tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping H2O_dimer_sep_pp tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

