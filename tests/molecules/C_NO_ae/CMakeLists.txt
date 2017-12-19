
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# C atom numerical orbital (SQD) tests, all electron
#

# Short tests
  LIST(APPEND CNOHFS_VMC_SCALARS "totenergy" "-37.68858615 0.002370") # total energy
  LIST(APPEND CNOSJS_VMC_SCALARS "totenergy" "-37.769812   0.001065") # total energy
  LIST(APPEND CNOSJS_VMC_SCALARS "variance"  "  0.632539   0.001819") # energy variance
  LIST(APPEND CNOSJS_DMC_SCALARS "totenergy" "-37.828225   0.000874") # total energy
  LIST(APPEND CNOSJS_DMC_SCALARS "variance"  "  0.622712   0.000844") # energy variance

  QMC_RUN_AND_CHECK(short-c_no-hf_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C_NO_ae"
                    C.NO_HF_short
                    C.NO_HF.short.in.xml
                    1 16
                    TRUE
                    0 CNOHFS_VMC_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-c_no-hf_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C_NO_ae"
                    C.NO_HF_short
                    C.NO_HF.short.in.xml
                    4 4
                    TRUE
                    0 CNOHFS_VMC_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-c_no-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C_NO_ae"
                    C.NO_J2_short
                    C.NO_J2.short.in.xml
                    1 16
                    TRUE
                    0 CNOSJS_VMC_SCALARS # VMC
                    2 CNOSJS_DMC_SCALARS # DMC
                    )

# Long tests
  LIST(APPEND CNOHFL_VMC_SCALARS "totenergy" "-37.68858615 0.000749") # total energy
  LIST(APPEND CNOSJL_VMC_SCALARS "totenergy" "-37.769812   0.000352") # total energy
  LIST(APPEND CNOSJL_VMC_SCALARS "variance"  "  0.632539   0.000600") # energy variance
  LIST(APPEND CNOSJL_DMC_SCALARS "totenergy" "-37.828225   0.000289") # total energy
  LIST(APPEND CNOSJL_DMC_SCALARS "variance"  "  0.622712   0.000279") # energy variance

  QMC_RUN_AND_CHECK(long-c_no-hf_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C_NO_ae"
                    C.NO_HF_long
                    C.NO_HF.long.in.xml
                    1 16
                    TRUE
                    0 CNOHFL_VMC_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-c_no-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C_NO_ae"
                    C.NO_J2_long
                    C.NO_J2.long.in.xml
                    1 16
                    TRUE
                    0 CNOSJL_VMC_SCALARS # VMC
                    2 CNOSJL_DMC_SCALARS # DMC
                    )

 ELSE()
    MESSAGE("Skipping C_NO_ae tests because the basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping C_NO_ae tests because the basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

