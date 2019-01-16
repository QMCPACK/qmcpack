
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# C atom numerical orbital (SQD) tests, all electron
#

# Short tests
  LIST(APPEND CNOHFS_VMC_SCALARS "totenergy" "-37.68858615 0.002286571486") # total energy
  LIST(APPEND CNOSJS_VMC_SCALARS "totenergy" "-37.769709615159   0.001025439060") # total energy
  LIST(APPEND CNOSJS_VMC_SCALARS "variance"  "  0.632569308194   0.002384473575") # energy variance
  LIST(APPEND CNOSJS_DMC_SCALARS "totenergy" "-37.828242297763   0.000839597533") # total energy
  LIST(APPEND CNOSJS_DMC_SCALARS "variance"  "  0.622617021499   0.000901570915") # energy variance

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
  LIST(APPEND CNOHFL_VMC_SCALARS "totenergy" "-37.68858615 0.000731146702") # total energy
  LIST(APPEND CNOSJL_VMC_SCALARS "totenergy" "-37.769709615159   0.000327891077") # total energy
  LIST(APPEND CNOSJL_VMC_SCALARS "variance"  "  0.632569308194   0.000762451558") # energy variance
  LIST(APPEND CNOSJL_DMC_SCALARS "totenergy" "-37.828242297763   0.000268138571") # total energy
  LIST(APPEND CNOSJL_DMC_SCALARS "variance"  "  0.622617021499   0.000287930738") # energy variance

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

