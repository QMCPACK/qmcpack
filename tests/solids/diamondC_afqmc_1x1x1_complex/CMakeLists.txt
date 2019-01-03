SET(TEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/choldump.h5)
SET(TEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/wfn.dat)

#symlink h5 file
MAYBE_SYMLINK(${TEST_HDF_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/choldump.h5)
#symlink wfn file
MAYBE_SYMLINK(${TEST_WFN_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/wfn.dat)

LIST(APPEND DIAMOND_AFQMC_SCALARS_CHOL "EnergyEstim__nume_real" "-379.2037 0.0282")

QMC_RUN_AND_CHECK("short-diamondC_afqmc_1x1x1_complex_cholesky"
                  "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
                  qmc_short_chol
                  qmc_short_chol.in.xml
                  16 1
                  TRUE
                  0 DIAMOND_AFQMC_SCALARS_CHOL
                  )

LIST(APPEND LONG_DIAMOND_AFQMC_SCALARS_CHOL "EnergyEstim__nume_real" "-379.2037 0.0043")

QMC_RUN_AND_CHECK("long-diamondC_afqmc_1x1x1_complex_cholesky"
                  "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
                  qmc_long_chol
                  qmc_long_chol.in.xml
                  16 1
                  TRUE
                  0 LONG_DIAMOND_AFQMC_SCALARS_CHOL
                  )

SET(TEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/thc_minimal/thc.h5)
SET(TEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/thc_minimal/wfn_thc.dat)

#symlink h5 file
MAYBE_SYMLINK(${TEST_HDF_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/thc.h5)
#symlink wfn file
MAYBE_SYMLINK(${TEST_WFN_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/wfn_thc.dat)

LIST(APPEND DIAMOND_AFQMC_SCALARS_THC "EnergyEstim__nume_real" "-10.33408 0.00432")

QMC_RUN_AND_CHECK("short-diamondC_afqmc_1x1x1_complex_thc"
                  "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
                  qmc_short_thc
                  qmc_short_thc.in.xml
                  16 1
                  TRUE
                  0 DIAMOND_AFQMC_SCALARS_THC
                  )

LIST(APPEND LONG_DIAMOND_AFQMC_SCALARS_THC "EnergyEstim__nume_real" "-10.33408 0.00143")

QMC_RUN_AND_CHECK("long-diamondC_afqmc_1x1x1_complex_thc"
                  "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
                  qmc_long_thc
                  qmc_long_thc.in.xml
                  16 1
                  TRUE
                  0 LONG_DIAMOND_AFQMC_SCALARS_THC
                  )
