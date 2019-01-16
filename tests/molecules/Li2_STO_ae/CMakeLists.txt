
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# Li2 molecular dimer STO tests, all electron
#   Reproduces UNR 1993
#

# Short tests
  LIST(APPEND LI2STOSJS_VMC_SCALARS "totenergy" "-14.947221704021  0.000491521212") # total energy
  LIST(APPEND LI2STOSJS_VMC_SCALARS "variance"  "  0.170589099144  0.001406847789") # energy variance
  LIST(APPEND LI2STOSJS_DMC_SCALARS "totenergy" "-14.989463714638  0.000690481867") # total energy
  LIST(APPEND LI2STOSJS_DMC_SCALARS "variance"  "  0.151150336595  0.000378095395") # energy variance

  QMC_RUN_AND_CHECK(short-li2_sto-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
                    Li2.STO.short
                    Li2.STO.short.in.xml
                    1 16
                    TRUE
                    0 LI2STOSJS_VMC_SCALARS # VMC
                    2 LI2STOSJS_DMC_SCALARS # DMC
                    )

  QMC_RUN_AND_CHECK(short-li2_sto-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
                    Li2.STO.short
                    Li2.STO.short.in.xml
                    4 4
                    TRUE
                    0 LI2STOSJS_VMC_SCALARS # VMC
                    2 LI2STOSJS_DMC_SCALARS # DMC
                    )

# Long tests
  LIST(APPEND LI2STOSJL_VMC_SCALARS "totenergy" "-14.947221704021  0.000157167233") # total energy
  LIST(APPEND LI2STOSJL_VMC_SCALARS "variance"  "  0.170589099144  0.000449849099") # energy variance
  LIST(APPEND LI2STOSJL_DMC_SCALARS "totenergy" "-14.989463714638  0.000220786248") # total energy
  LIST(APPEND LI2STOSJL_DMC_SCALARS "variance"  "  0.151150336595  0.000120898560") # energy variance

  QMC_RUN_AND_CHECK(long-li2_sto-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
                    Li2.STO.long
                    Li2.STO.long.in.xml
                    1 16
                    TRUE
                    0 LI2STOSJL_VMC_SCALARS # VMC
                    2 LI2STOSJL_DMC_SCALARS # DMC
                    )

 ELSE()
    MESSAGE("Skipping Li2_STO_ae tests because the basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping Li2_STO_ae tests because the basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

