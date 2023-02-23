#
# Homogeneous electron gas (HEG) tests
#

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
list(APPEND HEG14GSNI_SCALARS "totenergy" " 0.627711  0.000001") # total energy
list(APPEND HEG14GSNI_SCALARS "variance" " 0.000000  0.000001") # energy variance

qmc_run_and_check(
  short-heg_14_gamma-ni
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_NI
  heg-short-NI.xml
  1
  16
  TRUE
  0
  HEG14GSNI_SCALARS)

#
#  HEG14G - non-interacting Slater only DMC (exact)
#
list(APPEND HEG14GSNI_DMC_SCALARS "totenergy" " 0.627711  0.000001") # total energy
list(APPEND HEG14GSNI_DMC_SCALARS "variance" " 0.000000  0.000001") # energy variance

qmc_run_and_check(
  short-heg_14_gamma-ni_dmc
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_NI_dmc
  heg-short-NI-dmc.xml
  1
  16
  TRUE
  1
  HEG14GSNI_DMC_SCALARS # DMC
)

#
#  HEG14G - Slater only VMC (Hartree-Fock)
#
list(APPEND HEG14GSHF_SCALARS "totenergy" "-0.812484  0.000271") # total energy
list(APPEND HEG14GSHF_SCALARS "variance" " 0.193169  0.018612") # energy variance

qmc_run_and_check(
  short-heg_14_gamma-hf
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_HF
  heg-short-HF.xml
  1
  16
  TRUE
  0
  HEG14GSHF_SCALARS # VMC
)

#
#  HEG14G - Slater-Jastrow VMC
#
list(APPEND HEG14GSSJ_SCALARS "totenergy" "-1.073323  0.000271") # total energy
list(APPEND HEG14GSSJ_SCALARS "variance" " 0.024574  0.000181") # energy variance

qmc_run_and_check(
  short-heg_14_gamma-sj
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_SJ
  heg-short-SJ.xml
  1
  16
  TRUE
  0
  HEG14GSSJ_SCALARS # VMC
)

#
#  HEG14G - Slater-Jastrow VMC with newer input format
#
list(APPEND HEG14GSSJ_NEW_SCALARS "totenergy" "-1.073286  0.000271") # total energy
list(APPEND HEG14GSSJ_NEW_SCALARS "variance" " 0.024564  0.000342") # energy variance

# legacy driver
qmc_run_and_check(
  short-heg_14_gamma-sj_new
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_SJ_new
  heg-short-SJ-new.xml
  1
  16
  TRUE
  0
  HEG14GSSJ_NEW_SCALARS # VMC
)

# batch driver
qmc_run_and_check(
  short-heg_14_gamma-sj_batch
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma/batch"
  short_heg_SJ_batch
  short-heg-SJ-batch.xml
  1
  16
  TRUE
  0
  HEG14GSSJ_NEW_SCALARS # VMC
)

# batch deterministic
if(QMC_MIXED_PRECISION)
  list(APPEND DET_HEG14GSSJ_SCALARS "totenergy" "-1.24699388 0.00002")
  list(APPEND DET_HEG14GSSJ_SCALARS "kinetic" "0.77668380 0.000005")
  list(APPEND DET_HEG14GSSJ_SCALARS "potential" "-2.02367769 0.00002")
else()
  list(APPEND DET_HEG14GSSJ_SCALARS "totenergy" "-1.24699388 0.000001")
  list(APPEND DET_HEG14GSSJ_SCALARS "kinetic" "0.77668380 0.000001")
  list(APPEND DET_HEG14GSSJ_SCALARS "potential" "-2.02367769 0.000001")
endif()

qmc_run_and_check(
  deterministic-heg_14_gamma-sj-batch
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma/batch"
  deter_heg_SJ_batch
  deter-heg-SJ-batch.xml
  1
  1
  TRUE
  0
  DET_HEG14GSSJ_SCALARS # VMC
)

#
#  HEG14G - Slater-Jastrow DMC
#
list(APPEND HEG14GSSJ_DMC_SCALARS "totenergy" "-1.110199  0.000322") # total energy
list(APPEND HEG14GSSJ_DMC_SCALARS "variance" " 0.022894  0.000060") # energy variance

qmc_run_and_check(
  short-heg_14_gamma-sj_dmc
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_SJ_dmc
  heg-short-SJ-dmc.xml
  1
  16
  TRUE
  2
  HEG14GSSJ_DMC_SCALARS # DMC
)

#
# Long HEG 14 gamma tests
#

#
#  HEG14G - non-interacting Slater only VMC (exact)
#
list(APPEND HEG14GLNI_SCALARS "totenergy" " 0.627711  0.000001") # total energy
list(APPEND HEG14GLNI_SCALARS "variance" " 0.000000  0.000001") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-ni
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_NI
  heg-long-NI.xml
  1
  16
  TRUE
  0
  HEG14GLNI_SCALARS # VMC
)

#
#  HEG14G - non-interacting Slater only DMC (exact)
#
list(APPEND HEG14GLNI_DMC_SCALARS "totenergy" " 0.627711  0.000001") # total energy
list(APPEND HEG14GLNI_DMC_SCALARS "variance" " 0.000000  0.000001") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-ni_dmc
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_NI_dmc
  heg-long-NI-dmc.xml
  1
  16
  TRUE
  1
  HEG14GLNI_DMC_SCALARS # DMC
)

#
#  HEG14G - Slater only VMC (Hartree-Fock)
#
list(APPEND HEG14GLHF_SCALARS "totenergy" "-0.812484  0.000090") # total energy
list(APPEND HEG14GLHF_SCALARS "variance" " 0.193169  0.006142") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-hf
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_HF
  heg-long-HF.xml
  1
  16
  TRUE
  0
  HEG14GLHF_SCALARS # VMC
)

#
#  HEG14G - Slater-Jastrow VMC
#
list(APPEND HEG14GLSJ_SCALARS "totenergy" "-1.073323  0.000090") # total energy
list(APPEND HEG14GLSJ_SCALARS "variance" " 0.024574  0.000060") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-sj
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_SJ
  heg-long-SJ.xml
  1
  16
  TRUE
  0
  HEG14GLSJ_SCALARS # VMC
)

#
#  HEG14G - Slater-Jastrow VMC with newer input format
#
list(APPEND HEG14GLSJ_NEW_SCALARS "totenergy" "-1.073286  0.000090") # total energy
list(APPEND HEG14GLSJ_NEW_SCALARS "variance" " 0.024564  0.000113") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-sj_new
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_SJ_new
  heg-long-SJ-new.xml
  1
  16
  TRUE
  0
  HEG14GLSJ_NEW_SCALARS # VMC
)

#
#  HEG14G - Slater-Jastrow DMC
#
list(APPEND HEG14GLSJ_DMC_SCALARS "totenergy" "-1.110199  0.000106") # total energy
list(APPEND HEG14GLSJ_DMC_SCALARS "variance" " 0.022894  0.000020") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-sj_dmc
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_SJ_dmc
  heg-long-SJ-dmc.xml
  1
  16
  TRUE
  2
  HEG14GLSJ_DMC_SCALARS # DMC
)

#
#  HEG14G - Slater-Jastrow-Backflow VMC
#           Run 16x1, 4x4, 1x16 to test SJB cloning
#
list(APPEND HEG14GSSJB_SCALARS "totenergy" "-1.084963  0.0004") # total energy
list(APPEND HEG14GSSJB_SCALARS "variance" " 0.022667  0.000764") # energy variance

qmc_run_and_check(
  short-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_SJB
  heg-short-SJB.xml
  16
  1
  TRUE
  0
  HEG14GSSJB_SCALARS # VMC
)

qmc_run_and_check(
  short-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_SJB
  heg-short-SJB.xml
  4
  4
  TRUE
  0
  HEG14GSSJB_SCALARS # VMC
)

qmc_run_and_check(
  short-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_short_SJB
  heg-short-SJB.xml
  1
  16
  TRUE
  0
  HEG14GSSJB_SCALARS # VMC
)
#
#  HEG14G - Slater-Jastrow-Backflow VMC
#           Run 16x1, 4x4, 1x16 to test SJB cloning
#
list(APPEND HEG14GLSJB_SCALARS "totenergy" "-1.084963  0.000083") # total energy
list(APPEND HEG14GLSJB_SCALARS "variance" " 0.022667  0.000252") # energy variance

qmc_run_and_check(
  long-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_SJB
  heg-long-SJB.xml
  16
  1
  TRUE
  0
  HEG14GLSJB_SCALARS # VMC
)

qmc_run_and_check(
  long-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_SJB
  heg-long-SJB.xml
  4
  4
  TRUE
  0
  HEG14GLSJB_SCALARS # VMC
)

qmc_run_and_check(
  long-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  heg_long_SJB
  heg-long-SJB.xml
  1
  16
  TRUE
  0
  HEG14GLSJB_SCALARS # VMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_HEG14GSSJB_SCALARS "totenergy" "-0.92843602  0.00002")
  list(APPEND DET_HEG14GSSJB_SCALARS "kinetic" "0.76941696  0.000005")
  list(APPEND DET_HEG14GSSJB_SCALARS "potential" "-1.69785298 0.00002")
else()
  list(APPEND DET_HEG14GSSJB_SCALARS "totenergy" "-0.92845399  0.000001")
  list(APPEND DET_HEG14GSSJB_SCALARS "kinetic" "0.76941751  0.000001")
  list(APPEND DET_HEG14GSSJB_SCALARS "potential" "-1.69787150 0.000001")
endif()

qmc_run_and_check(
  deterministic-heg_14_gamma-sjb
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
  det_heg_short_SJB
  det_heg-short-SJB.xml
  1
  1
  TRUE
  0
  DET_HEG14GSSJB_SCALARS # VMC
)

if(NOT QMC_MIXED_PRECISION)
  simple_run_and_check(
    deterministic-heg_14_gamma-sjb-opt
    "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
    det_heg-short-SJB_opt.xml
    1
    1
    check_opt.py
    -s
    1
    -p
    det_heg_short_SJB
    -r
    ./qmc-ref/ref_det_heg_short_SJB.s001.opt.xml)

  list(APPEND DET_HEG14GSSJB_OPT_SCALARS "totenergy" "-0.93078813 0.000001")

  qmc_run_and_check(
    deterministic-heg_14_gamma-sjb-opt_vmc
    "${qmcpack_SOURCE_DIR}/tests/heg/heg_14_gamma"
    det_heg_short_SJB
    det_heg-short-SJB_opt.xml
    1
    1
    TRUE
    2
    DET_HEG14GSSJB_OPT_SCALARS # VMC
  )

endif()
