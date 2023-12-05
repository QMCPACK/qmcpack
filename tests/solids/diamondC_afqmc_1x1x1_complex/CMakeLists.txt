# RHF CHOLESKY
set(TEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/ham_chol_sc.h5)
set(TEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/wfn_rhf.dat)

#symlink h5 file
maybe_symlink(${TEST_HDF_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/choldump.h5)
#symlink wfn file
maybe_symlink(${TEST_WFN_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/wfn.dat)

list(APPEND DIAMOND_AFQMC_SCALARS_CHOL "EnergyEstim__nume_real" "-10.55201 0.00332")

qmc_run_and_check(
  "short-diamondC_afqmc_1x1x1_complex_cholesky"
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
  qmc_short_chol
  qmc_short_chol.in.xml
  16
  1
  TRUE
  0
  DIAMOND_AFQMC_SCALARS_CHOL)

list(APPEND LONG_DIAMOND_AFQMC_SCALARS_CHOL "EnergyEstim__nume_real" "-10.55201 0.00109")

qmc_run_and_check(
  "long-diamondC_afqmc_1x1x1_complex_cholesky"
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
  qmc_long_chol
  qmc_long_chol.in.xml
  16
  1
  TRUE
  0
  LONG_DIAMOND_AFQMC_SCALARS_CHOL)

# RHF ISDF/THC
set(TEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_szv/ham_thc_sc.h5)
set(TEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_szv/wfn_rhf.dat)

#symlink h5 file
maybe_symlink(${TEST_HDF_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/thc.h5)
#symlink wfn file
maybe_symlink(${TEST_WFN_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/wfn_thc.dat)

list(APPEND DIAMOND_AFQMC_SCALARS_THC "EnergyEstim__nume_real" "-10.33408 0.00432")

qmc_run_and_check(
  "short-diamondC_afqmc_1x1x1_complex_thc"
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
  qmc_short_thc
  qmc_short_thc.in.xml
  16
  1
  TRUE
  0
  DIAMOND_AFQMC_SCALARS_THC)

list(APPEND LONG_DIAMOND_AFQMC_SCALARS_THC "EnergyEstim__nume_real" "-10.33408 0.00143")

qmc_run_and_check(
  "long-diamondC_afqmc_1x1x1_complex_thc"
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
  qmc_long_thc
  qmc_long_thc.in.xml
  16
  1
  TRUE
  0
  LONG_DIAMOND_AFQMC_SCALARS_THC)

# UHF CHOLESKY
# Need to remake benchmark following change in hamiltonian/wavefunction.
#SET(TEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_uhf/choldump.h5)
#SET(TEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_uhf/wfn.dat)

##symlink h5 file
#MAYBE_SYMLINK(${TEST_HDF_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/choldump_uhf.h5)
##symlink wfn file
#MAYBE_SYMLINK(${TEST_WFN_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/wfn_uhf.dat)

#LIST(APPEND DIAMOND_AFQMC_SCALARS_UHF "EnergyEstim__nume_real" "-10.55151 0.00312")

#QMC_RUN_AND_CHECK("short-diamondC_afqmc_1x1x1_complex_uhf"
#"${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
#qmc_short_uhf
#qmc_short_uhf.in.xml
#16 1
#TRUE
#0 DIAMOND_AFQMC_SCALARS_UHF
#)

#LIST(APPEND LONG_DIAMOND_AFQMC_SCALARS_UHF "EnergyEstim__nume_real" "-10.55151 0.00103")

#QMC_RUN_AND_CHECK("long-diamondC_afqmc_1x1x1_complex_uhf"
#"${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
#qmc_long_uhf
#qmc_long_uhf.in.xml
#16 1
#TRUE
#0 LONG_DIAMOND_AFQMC_SCALARS_UHF
#)

# PHMSD CHOLESKY
#SET(TEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/choldump_phmsd.h5)
#SET(TEST_WFN_INPUT ${qmcpack_SOURCE_DIR}/tests/afqmc/C_rhf/wfn_phmsd.dat)
##symlink h5 file
#MAYBE_SYMLINK(${TEST_HDF_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/choldump_phmsd.h5)
##symlink wfn file
#MAYBE_SYMLINK(${TEST_WFN_INPUT} ${CMAKE_CURRENT_BINARY_DIR}/wfn_phmsd.dat)
#
#LIST(APPEND DIAMOND_AFQMC_SCALARS_PHMSD "EnergyEstim__nume_real" "-10.55085 0.00368")
#
#QMC_RUN_AND_CHECK("short-diamondC_afqmc_1x1x1_complex_phmsd"
#                  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
#                  qmc_short_phmsd
#                  qmc_short_phmsd.in.xml
#                  16 1
#                  TRUE
#                  0 DIAMOND_AFQMC_SCALARS_PHMSD
#                  )
#
#LIST(APPEND LONG_DIAMOND_AFQMC_SCALARS_PHMSD "EnergyEstim__nume_real" "-10.55085 0.00119")
#
#QMC_RUN_AND_CHECK("long-diamondC_afqmc_1x1x1_complex_phmsd"
#                  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_afqmc_1x1x1_complex"
#                  qmc_long_phmsd
#                  qmc_long_phmsd.in.xml
#                  16 1
#                  TRUE
#                  0 LONG_DIAMOND_AFQMC_SCALARS_PHMSD
#                  )
