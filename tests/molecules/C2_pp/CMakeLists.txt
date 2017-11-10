
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

# C2 TZ PP multi determinants
# Short tests
  LIST(APPEND C2MSDJS_VMC_SCALARS "totenergy" "-11.064964  0.000865") # total energy
  LIST(APPEND C2MSDJS_VMC_SCALARS "variance"  "  0.149682  0.003928") # energy variance

  QMC_RUN_AND_CHECK(short-C2_pp-msdj_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_pp"
                    C2_TZ_PP_MSDJ
                    qmc_short_pp_msdj.xml
                    1 16
                    ${MP_SUCCESS}
                    1 C2MSDJS_VMC_SCALARS # VMC
                    )

  LIST(APPEND C2MSDJL_VMC_SCALARS "totenergy" "-11.064964  0.000274") # total energy
  LIST(APPEND C2MSDJL_VMC_SCALARS "variance"  "  0.149682  0.001242") # energy variance

  QMC_RUN_AND_CHECK(long-C2_pp-msdj_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_pp"
                    C2_TZ_PP_MSDJ
                    qmc_long_pp_msdj.xml
                    1 16
                    ${MP_SUCCESS}
                    1 C2MSDJL_VMC_SCALARS # VMC
                    )
 ELSE()
    MESSAGE("Skipping C2_pp tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping C2_pp tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

