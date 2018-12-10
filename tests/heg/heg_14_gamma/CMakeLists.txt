

#
# Homogeneous electron gas (HEG) tests
#
IF (NOT QMC_CUDA)

#
# Tests for HEG with 14 electrons at the gamma point and various wavefunctions
# See tests/heg/heg_14_gamma/README for full details
#

#
# Short HEG 14 gamma tests
#

#
#  HEG14G - non-interacting Slater only VMC (exact)
#
  LIST(APPEND HEG14GSNI_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GSNI_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-ni
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_NI
                    heg-short-NI.xml
                    1 16
                    TRUE
                    0 HEG14GSNI_SCALARS
                    )

#
#  HEG14G - non-interacting Slater only DMC (exact)
#
  LIST(APPEND HEG14GSNI_DMC_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GSNI_DMC_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-ni_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_NI_dmc
                    heg-short-NI-dmc.xml
                    1 16
                    TRUE
                    1 HEG14GSNI_DMC_SCALARS # DMC
                    )

#
#  HEG14G - Slater only VMC (Hartree-Fock)
#
  LIST(APPEND HEG14GSHF_SCALARS "totenergy" "-0.812484  0.000271") # total energy
  LIST(APPEND HEG14GSHF_SCALARS "variance"  " 0.193169  0.018612") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-hf
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_HF
                    heg-short-HF.xml
                    1 16
                    TRUE
                    0 HEG14GSHF_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC
#
  LIST(APPEND HEG14GSSJ_SCALARS "totenergy" "-1.073323  0.000271") # total energy
  LIST(APPEND HEG14GSSJ_SCALARS "variance"  " 0.024574  0.000181") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sj
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJ
                    heg-short-SJ.xml
                    1 16
                    TRUE
                    0 HEG14GSSJ_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC with newer input format
#
  LIST(APPEND HEG14GSSJ_NEW_SCALARS "totenergy" "-1.073286  0.000271") # total energy
  LIST(APPEND HEG14GSSJ_NEW_SCALARS "variance"  " 0.024564  0.000342") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sj_new
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJ_new
                    heg-short-SJ-new.xml
                    1 16
                    TRUE
                    0 HEG14GSSJ_NEW_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow DMC
#
  LIST(APPEND HEG14GSSJ_DMC_SCALARS "totenergy" "-1.110199  0.000322") # total energy
  LIST(APPEND HEG14GSSJ_DMC_SCALARS "variance"  " 0.022894  0.000060") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJ_dmc
                    heg-short-SJ-dmc.xml
                    1 16
                    TRUE
                    2 HEG14GSSJ_DMC_SCALARS # DMC
                    )

#
# Long HEG 14 gamma tests
#

#
#  HEG14G - non-interacting Slater only VMC (exact)
#
  LIST(APPEND HEG14GLNI_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GLNI_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-ni
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_NI
                    heg-long-NI.xml
                    1 16
                    TRUE
                    0 HEG14GLNI_SCALARS # VMC
                    )

#
#  HEG14G - non-interacting Slater only DMC (exact)
#
  LIST(APPEND HEG14GLNI_DMC_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GLNI_DMC_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-ni_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_NI_dmc
                    heg-long-NI-dmc.xml
                    1 16
                    TRUE
                    1 HEG14GLNI_DMC_SCALARS # DMC
                    )

#
#  HEG14G - Slater only VMC (Hartree-Fock)
#
  LIST(APPEND HEG14GLHF_SCALARS "totenergy" "-0.812484  0.000090") # total energy
  LIST(APPEND HEG14GLHF_SCALARS "variance"  " 0.193169  0.006142") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-hf
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_HF
                    heg-long-HF.xml
                    1 16
                    TRUE
                    0 HEG14GLHF_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC
#
  LIST(APPEND HEG14GLSJ_SCALARS "totenergy" "-1.073323  0.000090") # total energy
  LIST(APPEND HEG14GLSJ_SCALARS "variance"  " 0.024574  0.000060") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sj
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJ
                    heg-long-SJ.xml
                    1 16
                    TRUE
                    0 HEG14GLSJ_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC with newer input format
#
  LIST(APPEND HEG14GLSJ_NEW_SCALARS "totenergy" "-1.073286  0.000090") # total energy
  LIST(APPEND HEG14GLSJ_NEW_SCALARS "variance"  " 0.024564  0.000113") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sj_new
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJ_new
                    heg-long-SJ-new.xml
                    1 16
                    TRUE
                    0 HEG14GLSJ_NEW_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow DMC
#
  LIST(APPEND HEG14GLSJ_DMC_SCALARS "totenergy" "-1.110199  0.000106") # total energy
  LIST(APPEND HEG14GLSJ_DMC_SCALARS "variance"  " 0.022894  0.000020") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJ_dmc
                    heg-long-SJ-dmc.xml
                    1 16
                    TRUE
                    2 HEG14GLSJ_DMC_SCALARS # DMC
                    )

IF (NOT QMC_COMPLEX)
#
#  HEG14G - Slater-Jastrow-Backflow VMC
#           Run 16x1, 4x4, 1x16 to test SJB cloning
#
  LIST(APPEND HEG14GSSJB_SCALARS "totenergy" "-1.084963  0.000251") # total energy
  LIST(APPEND HEG14GSSJB_SCALARS "variance"  " 0.022667  0.000764") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJB
                    heg-short-SJB.xml
                    16 1
                    TRUE
                    0 HEG14GSSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJB
                    heg-short-SJB.xml
                    4 4
                    TRUE
                    0 HEG14GSSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJB
                    heg-short-SJB.xml
                    1 16
                    TRUE
                    0 HEG14GSSJB_SCALARS # VMC
                    )
#
#  HEG14G - Slater-Jastrow-Backflow VMC
#           Run 16x1, 4x4, 1x16 to test SJB cloning
#
  LIST(APPEND HEG14GLSJB_SCALARS "totenergy" "-1.084963  0.000083") # total energy
  LIST(APPEND HEG14GLSJB_SCALARS "variance"  " 0.022667  0.000252") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJB
                    heg-long-SJB.xml
                    16 1
                    TRUE
                    0 HEG14GLSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJB
                    heg-long-SJB.xml
                    4 4
                    TRUE
                    0 HEG14GLSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJB
                    heg-long-SJB.xml
                    1 16
                    TRUE
                    0 HEG14GLSJB_SCALARS # VMC
                    )

 ELSE()
    MESSAGE("Skipping HEG backflow tests because they are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()

ELSE()
    MESSAGE("Skipping all HEG tests because they are not supported by CUDA build (QMC_CUDA=1)")  
ENDIF()

