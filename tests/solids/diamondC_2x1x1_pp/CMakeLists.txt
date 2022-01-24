list(APPEND DIAMOND2_SCALARS "totenergy" "-21.668293 0.007")
list(APPEND DIAMOND2_SCALARS "kinetic" "20.065151 0.044")
list(APPEND DIAMOND2_SCALARS "potential" "-41.733377 0.045")
list(APPEND DIAMOND2_SCALARS "eeenergy" "-5.428561 0.012")
list(APPEND DIAMOND2_SCALARS "ionion" "-25.5513269 0.001")
list(APPEND DIAMOND2_SCALARS "localecp" "-12.29402 0.051")
list(APPEND DIAMOND2_SCALARS "nonlocalecp" "1.540530 0.018")
list(APPEND DIAMOND2_SCALARS "samples" "32000 0.0")
#  LIST(APPEND DIAMOND2_SCALARS "flux" "0.0 0.7")

qmc_run_and_check(
  short-diamondC_2x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  DIAMOND2_SCALARS # VMC
)
if(NOT QMC_CUDA)
  # 16 batches of 1 walker each
  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short
    qmc_short_vmcbatch.in.xml
    1
    16
    TRUE
    0
    DIAMOND2_SCALARS # VMC
  )

  # 4 batches of 4 walkers
  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch-4b_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short
    qmc_short_vmcbatch_4b.in.xml
    1
    4
    TRUE
    0
    DIAMOND2_SCALARS # VMC
  )

  # 4 batches of 4 walkers with asynchronous WF computation
  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch-4b-asynctwf_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short
    qmc_short_vmcbatch_4b_asynctwf.in.xml
    1
    4
    TRUE
    0
    DIAMOND2_SCALARS # VMC
  )

  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch-4b_sdbatch_jas
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short
    qmc_short_sdbatch_vmcbatch_4b.in.xml
    1
    4
    TRUE
    0
    DIAMOND2_SCALARS # VMC
  )
endif()
qmc_run_and_check(
  short-diamondC_2x1x1_pp-delayed_update-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short
  qmc_short_delayedupdate.in.xml
  1
  16
  TRUE
  0
  DIAMOND2_SCALARS # VMC
)

qmc_run_and_check(
  short-diamondC_2x1x1_hybridrep_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short
  qmc_short_hybridrep.in.xml
  1
  16
  ${SUCCESS_STATUS_CUDA}
  0
  DIAMOND2_SCALARS # VMC
)

# Reference DMC run in qmc-ref "-21.844975 0.001113"
list(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
qmc_run_and_check(
  short-diamondC_2x1x1_pp-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short_vmc_dmc
  qmc_short_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  DIAMOND2_DMC_SCALARS # DMC
)

# 1 rank 4 threads 16 walkers (4 crowds of 4 walkers)
list(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
qmc_run_and_check(
  short-diamondC_2x1x1_pp-dmc_4r_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short_vmc_dmc_4r
  qmc_short_vmc_dmc_4r.in.xml
  1
  4
  TRUE
  1
  DIAMOND2_DMC_SCALARS # DMC
)

if(NOT QMC_CUDA)
  #16 batches of 16 walkers in one rank
  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch_dmcbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short_vmcbatch_dmcbatch
    qmc_short_vmcbatch_dmcbatch.in.xml
    1
    16
    TRUE
    1
    DIAMOND2_DMC_SCALARS # DMC
  )

  #4 batches of 16 walkers in one rank, batched slater det
  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch_dmcbatch-4b_sdbatch_jas
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short_sdbatch_vmcbatch_dmcbatch
    qmc_short_sdbatch_vmcbatch_dmcbatch.in.xml
    1
    4
    TRUE
    1
    DIAMOND2_DMC_SCALARS # DMC
  )

  # 4 batches of 16 walkers over 4 ranks
  qmc_run_and_check(
    short-diamondC_2x1x1_pp-vmcbatch_dmcbatch_4r_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    qmc_short_vmcbatch_dmcbatch_4r
    qmc_short_vmcbatch_dmcbatch_4r.in.xml
    4
    4
    TRUE
    1
    DIAMOND2_DMC_SCALARS # DMC
  )
endif()

list(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
qmc_run_and_check(
  short-diamondC_2x1x1_pp-dmc-reconf_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short_vmc_dmc
  qmc_short_vmc_dmc_reconf.in.xml
  1
  16
  TRUE
  1
  DIAMOND2_DMC_SCALARS # DMC
)

# Reference DMC run as qmc-ref -21.858821 +/- 0.001207 with 2500 blocks
list(APPEND DIAMOND2_DMC_DLA_SCALARS "totenergy" "-21.858821 0.012")
qmc_run_and_check(
  short-diamondC_2x1x1_pp-dmc-dla_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short_vmc_dmc_dla
  qmc_short_vmc_dmc_dla.in.xml
  1
  16
  TRUE
  1
  DIAMOND2_DMC_DLA_SCALARS # DMC
)

# Excited state

list(APPEND DIAMOND2_EXCITED_SCALARS "totenergy" "-21.117471 0.007937")
list(APPEND DIAMOND2_EXCITED_SCALARS "kinetic" "20.070595 0.043211")
list(APPEND DIAMOND2_EXCITED_SCALARS "potential" "-41.188066 0.044583")
list(APPEND DIAMOND2_EXCITED_SCALARS "eeenergy" "-5.453806 0.013619")
list(APPEND DIAMOND2_EXCITED_SCALARS "ionion" "-25.551327 0.00001")
list(APPEND DIAMOND2_EXCITED_SCALARS "localecp" "-11.732568 0.050756")
list(APPEND DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.549634 0.018029")
list(APPEND DIAMOND2_EXCITED_SCALARS "samples" "32000 0.0")

qmc_run_and_check(
  short-diamondC_2x1x1_pp-vmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short_excited
  qmc_short_excited.in.xml
  1
  16
  TRUE
  0
  DIAMOND2_EXCITED_SCALARS # VMC
)

list(APPEND DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-21.359004 0.012291")

qmc_run_and_check(
  short-diamondC_2x1x1_pp-dmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_short_vmc_dmc_excited
  qmc_short_vmc_dmc_excited.in.xml
  1
  16
  TRUE
  1
  DIAMOND2_DMC_EXCITED_SCALARS # DMC
)

#
# Long tests
#

list(APPEND LONG_DIAMOND2_SCALARS "totenergy" "-21.668293 0.00023")
list(APPEND LONG_DIAMOND2_SCALARS "samples" "30720000 0.0")
list(APPEND LONG_DIAMOND2_SCALARS "flux" "0.0 0.05")

qmc_run_and_check(
  long-diamondC_2x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND2_SCALARS # VMC
)

# Reference DMC run in qmc-ref "-21.844975 0.001113"
list(APPEND LONG_DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.002879")
qmc_run_and_check(
  long-diamondC_2x1x1_pp-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_long_vmc_dmc
  qmc_long_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND2_DMC_SCALARS # DMC
)

qmc_run_and_check(
  long-diamondC_2x1x1_pp-dmc-reconf_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_long_vmc_dmc
  qmc_long_vmc_dmc_reconf.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND2_DMC_SCALARS # DMC
)

# Reference DMC run as qmc-ref -21.858821 +/- 0.001207 with 2500 blocks
list(APPEND LONG_DIAMOND2_DMC_DLA_SCALARS "totenergy" "-21.858821 0.0038")
qmc_run_and_check(
  long-diamondC_2x1x1_pp-dmc-dla_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_long_vmc_dmc_dla
  qmc_long_vmc_dmc_dla.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND2_DMC_DLA_SCALARS # DMC
)

# Excited state
list(APPEND LONG_DIAMOND2_EXCITED_SCALARS "totenergy" "-21.117471 0.00027")
list(APPEND LONG_DIAMOND2_EXCITED_SCALARS "samples" "30720000 0.0")
list(APPEND LONG_DIAMOND2_EXCITED_SCALARS "flux" "-0.042346 0.075078")

qmc_run_and_check(
  long-diamondC_2x1x1_pp-vmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_long_excited
  qmc_long_excited.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND2_EXCITED_SCALARS # VMC
)

list(APPEND LONG_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-21.359004 0.004056")
qmc_run_and_check(
  long-diamondC_2x1x1_pp-dmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
  qmc_long_vmc_dmc_excited
  qmc_long_vmc_dmc_excited.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND2_DMC_EXCITED_SCALARS # DMC
)

#
# Deterministic tests
#

# VMC RUN
if(NOT QMC_CUDA)
  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-22.30282849 0.00003531")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.25970008 0.00012922")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.56252857 0.00009881")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-7.25690848 0.00002388")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55133296 0.0000014")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-10.85086697 0.00021079")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.09657984 0.00008703")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "kecorr" "0.030264591387 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-5.34650349 0.00047629")
    else()
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-22.29892024 0.00008588")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.26598410 0.00014541")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.56490434 0.00007327")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-7.25613581 0.00002740")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55133291 0.00000135")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-10.85690056 0.00018856")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.09946494 0.00010037")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "kecorr" "0.030264591387 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-5.36613862 0.00026178")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-22.30283274 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.25965947 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.56249222 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-7.25691718 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55132676 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-10.85078393 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.09653564 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "kecorr" "0.030264591387 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-5.34634371 0.000001")
    else()
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-22.29890271 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.26597773 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.56488044 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-7.25613869 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55132676 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-10.85686146 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.09944647 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "kecorr" "0.030264591387 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-5.36609515 0.000001")
    endif()

    #  QMC_RUN_AND_CHECK(deterministic-diamondC_2x1x1_pp-vmcbatch_sdj
    #                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    #                    det_qmc_short
    #                    det_qmc_short_vmcbatch.in.xml
    #                    1 1
    #                    TRUE
    #                    0 DET_DIAMOND2_SCALARS # VMC
    #                    )

  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short
    det_qmc_short.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_SCALARS # VMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-21.82533046 0.00003607")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "22.21682064 0.00005723")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-44.04215109 0.00004795")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-6.22921220 0.00001674")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55133296 0.0000014")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-13.30257077 0.00005318")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.04096483 0.00000726")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-19.91569882 0.00047829")
    else()
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-21.82258864 0.00002587")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "22.21966036 0.00006011")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-44.04224900 0.00004784")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-6.22922160 0.00001946")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55133291 0.00000135")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-13.30271625 0.00004867")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.04102175 0.00000802")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-19.92297713 0.0005331")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-21.82533699 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "22.21682105 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-44.04215804 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-6.22920728 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55132710 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-13.30258890 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.04096524 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-19.91551775 0.000001")
    else()
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-21.82258916 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "22.21966742 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-44.04225658 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-6.22921778 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55132710 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-13.30273399 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.04102228 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-19.92285549 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmc_sdj_excited
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_excited
    det_qmc_short_excited.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_EXCITED_SCALARS # VMC
  )

  # VMC DELAYED-UPDATE
  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-20.94751794 0.00001857")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "16.53206588 0.00002208")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-37.47958111 0.00002475")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-4.25515637 0.00000418")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55133296 0.0000014")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-9.39410011 0.00002983")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "1.72100833 0.00000237")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-12.47276797 0.00007329")
    else()
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-20.94387477 0.0000242")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "16.53548278 0.00002493")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-37.47935755 0.00002767")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-4.25504497 0.00000506")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55133291 0.00000135")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-9.39386018 0.00002929")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "1.72088050 0.0000045")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-12.48033341 0.00006604")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-20.94752145 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "16.53205375 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-37.47957520 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-4.25515766 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55132676 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-9.39410085 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "1.72101007 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-12.47274689 0.000001")
    else()
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-20.94387935 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "16.53547309 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-37.47935243 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-4.25504502 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55132676 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-9.39386209 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "1.72088142 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-12.48031692 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-delayed_update-vmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short
    det_qmc_short_delayedupdate.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_DELAY_SCALARS # VMC
  )

  # DMC RUN
  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-21.04525299 0.00001029")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "22.67695498 0.00001003")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-43.72220797 0.00000515")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-6.17624904 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55133296 0.0000014")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-9.10928079 0.00000442")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "-2.88534519 0.00000177")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-21.05326001 0.00001462")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "22.66898761 0.00001749")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-43.72224762 0.00000912")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-6.17623877 0.00000144")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55133291 0.00000135")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-9.10933337 0.00000562")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "-2.88534257 0.00000442")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-21.04525132 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "22.67697372 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-43.72222505 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-6.17624438 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55132676 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-9.10929967 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "-2.88535424 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-21.05325436 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "22.66901302 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-43.72226738 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-6.17624150 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55132676 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-9.10934270 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "-2.88535642 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-dmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmc_dmc
    det_qmc_short_vmc_dmc.in.xml
    1
    1
    TRUE
    1
    DET_DIAMOND2_DMC_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-19.57990811 0.00002128")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "24.19436676 0.00002563")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.77427487 0.00000861")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-5.85898725 0.00000192")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55133296 0.0000014")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.64208323 0.00000777")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "4.27812857 0.00001229")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-19.58550132 0.00001859")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "24.18875877 0.00002343")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.77426009 0.00001142")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-5.85898179 0.00000163")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55133291 0.00000135")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.64210613 0.00001074")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "4.27816074 0.00001354")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-19.57992731 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "24.19434492 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.77427223 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-5.85898752 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55132710 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.64208258 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "4.27812496 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-19.58550609 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "24.18875130 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.77425738 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-5.85898144 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55132710 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.64210491 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "4.27815606 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-dmc_sdj_excited
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmc_dmc_excited
    det_qmc_short_vmc_dmc_excited.in.xml
    1
    1
    TRUE
    1
    DET_DIAMOND2_DMC_EXCITED_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      # Cycle 1
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "totenergy" "-22.05192916 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "kinetic" "19.21103034 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "potential" "-41.26295939 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "eeenergy" "-6.68405093 0.00003")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "localecp" "-10.92577863 0.0002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "nonlocalecp" "1.89820362 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "flux" "-0.31336816 0.001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "samples" "18.00000000 0.0")
      #Cycle 2
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "totenergy" "-21.30130291 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "kinetic" "19.77013829 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "potential" "-41.07144136 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "eeenergy" "-6.35883147 0.00003")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "localecp" "-10.53491810 0.0002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "nonlocalecp" "1.37364166 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "flux" "-9.80945206 0.001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "samples" "18.00000000 0.0")
    else()
      #Cycle 1
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "totenergy" "-22.05264378 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "kinetic" "19.21031605 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "potential" "-41.26295995 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "eeenergy" "-6.68405093 0.00003")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "localecp" "-10.92577863 0.0002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "nonlocalecp" "1.89820306 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "flux" "-0.31142934 0.001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "samples" "18.00000000 0.0")
      #Cycle 2
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "totenergy" "-21.31027285 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "kinetic" "19.76116869 0.0005")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "potential" "-41.07144191 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "eeenergy" "-6.35883147 0.00003")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "localecp" "-10.53491810 0.0002")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "nonlocalecp" "1.37364111 0.0001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "flux" "-9.79154968 0.001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "samples" "18.00000000 0.0")
    endif()
  else()
    if(QMC_COMPLEX)
      #Cycle 1
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "totenergy" "-22.05191472 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "kinetic" "19.21103957 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "potential" "-41.26295429 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "eeenergy" "-6.68404915 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "localecp" "-10.92578832 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "nonlocalecp" "1.89821036 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "flux" "-0.31336458 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "samples" "18.00000000 0.0")
      #Cycle 2
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "totenergy" "-21.30129940 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "kinetic" "19.77014450 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "potential" "-41.07144390 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "eeenergy" "-6.35882865 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "localecp" "-10.53492678 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "nonlocalecp" "1.37363872 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "flux" "-9.80949941 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "samples" "18.00000000 0.0")
    else()
      #Cycle 1
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "totenergy" "-22.05263010 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "kinetic" "19.21032487 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "potential" "-41.26295497 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "eeenergy" "-6.68404915 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "localecp" "-10.92578832 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "nonlocalecp" "1.89820968 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "flux" "-0.31141016 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS "samples" "18.00000000 0.0")
      #Cycle 2
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "totenergy" "-21.31026809 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "kinetic" "19.76117659 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "potential" "-41.07144467 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "eeenergy" "-6.35882865 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "localecp" "-10.53492678 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "nonlocalecp" "1.37363795 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "flux" "-9.79159774 0.000001")
      list(APPEND DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS "samples" "18.00000000 0.0")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_mwalkers
    det_qmc_short_vmcbatch_mwalkers.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS # VMC cycle 1
    1
    DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS # VMC cycle 2
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-mwalkers_sdbatch_jas
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_mwalkers.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_VMC_BATCH_MWALKERS1_SCALARS # VMC cycle 1
    1
    DET_DIAMOND2_VMC_BATCH_MWALKERS2_SCALARS # VMC cycle 2
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy" "-21.75695352 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic" "21.58412194 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential" "-43.34107537 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy" "-4.93051914 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp" "-14.55016309 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp" "1.69094030 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "flux" "14.56015197 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy" "-21.40819399 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic" "21.75835414 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential" "-43.16654810 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy" "-5.31127989 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp" "-13.83088578 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp" "1.52695101 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "flux" "1.55142989 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy" "-21.76355918 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic" "24.31653559 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential" "-46.08009462 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy" "-5.15140370 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp" "-17.32746386 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp" "1.95010638 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "flux" "-15.81523996 0.001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy" "-21.76487128 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic" "20.67891363 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential" "-42.44378484 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy" "-4.97687774 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp" "-13.97216887 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp" "2.05659522 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "flux" "22.16100872 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy" "-21.34835919 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic" "22.16238919 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential" "-43.51074829 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy" "-5.07821053 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp" "-14.40620204 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp" "1.52499772 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "flux" "1.35192553 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy" "-21.72500894 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic" "24.02746675 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential" "-45.75247577 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy" "-5.05948524 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp" "-16.92160668 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp" "1.77994960 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "flux" "-15.34938772 0.001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy" "-21.77970598 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic" "22.07862416 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential" "-43.85833014 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy" "-4.98351398 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp" "-14.61422806 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp" "1.29073910 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "flux" "-10.01746831 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy" "-21.64880981 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic" "21.77217213 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential" "-43.42098194 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy" "-4.73027440 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp" "-14.29275394 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp" "1.15337360 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "flux" "-1.54466089 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy" "-21.45587606 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic" "20.81022717 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential" "-42.26610323 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy" "-4.20981722 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp" "-12.61049505 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp" "0.10553624 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "flux" "-4.19778396 0.000001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy" "-21.77684587 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic" "22.08147768 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential" "-43.85832356 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy" "-4.98351398 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp" "-14.61422806 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp" "1.29074568 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "flux" "-10.02287436 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy" "-21.64788214 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic" "21.77346471 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential" "-43.42134685 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy" "-4.73021029 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp" "-14.29312925 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp" "1.15331988 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "flux" "-1.54298684 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy" "-21.45991354 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic" "20.80700067 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential" "-42.26691421 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy" "-4.20986831 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp" "-12.61110393 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp" "0.10538522 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "flux" "-4.18809634 0.00001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_vmcbatch_dmcbatch_mwalkers.in.xml
    1
    4
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers-serialized_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_mwalkers_serialized
    det_qmc_short_vmcbatch_dmcbatch_mwalkers_serialized.in.xml
    1
    4
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers.in.xml
    1
    4
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers-checks_sdbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers_with_checks.in.xml
    1
    4
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy" "-21.58341272 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic" "19.37229488 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential" "-40.95570760 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy" "-5.29734575 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp" "-11.85045674 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp" "1.74342834 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "flux" "-6.02205066 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy" "-21.76467747 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic" "19.49291275 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential" "-41.25759022 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy" "-5.74480542 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp" "-11.45744019 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp" "1.49598884 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "flux" "-1.28278001 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy" "-22.22433395 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic" "19.10457607 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential" "-41.32891006 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy" "-5.53681243 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp" "-12.48670659 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp" "2.24594241 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "flux" "-6.00048904 0.001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy" "-21.58530618 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic" "19.37039175 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential" "-40.95569793 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy" "-5.29734575 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp" "-11.85045674 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp" "1.74343801 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "flux" "-6.01896369 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy" "-21.76315171 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic" "19.49451168 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential" "-41.25766339 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy" "-5.74484542 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp" "-11.45746322 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp" "1.49597869 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "flux" "-1.29073155 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy" "-22.22611131 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic" "19.10316757 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential" "-41.32927901 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy" "-5.53677170 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion" "-25.55133345 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp" "-12.48674850 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp" "2.24557464 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "flux" "-5.99726842 0.001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy" "-21.69113073 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic" "19.62048294 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential" "-41.31161367 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy" "-5.78662718 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp" "-12.25902950 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp" "2.28537020 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "flux" "-2.60864883 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy" "-21.76398145 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic" "20.42744622 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential" "-42.19142767 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy" "-5.74796574 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp" "-12.53877422 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp" "1.64663949 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "flux" "-6.56613667 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy" "-21.96323690 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic" "19.86538294 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential" "-41.82861984 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy" "-5.92974250 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp" "-12.66622927 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp" "2.31867912 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "flux" "-5.04750235 0.000001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy" "-21.69329544 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic" "19.61830510 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential" "-41.31160054 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy" "-5.78662718 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp" "-12.25902950 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp" "2.28538333 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "flux" "-2.60446436 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples" "144.0000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy" "-21.76448475 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic" "20.42746547 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential" "-42.19195022 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy" "-5.74793018 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp" "-12.53946292 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp" "1.64677007 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "flux" "-6.56446708 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy" "-21.96321405 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic" "19.86642474 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential" "-41.82963879 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy" "-5.92973716 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp" "-12.66715052 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp" "2.31857608 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "flux" "-5.04914162 0.00001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_vmcbatch_dmcbatch_mwalkers.in.xml
    1
    16
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers.in.xml
    1
    16
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy" "-21.58019753 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic" "19.06962113 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential" "-40.64981855 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy" "-4.95910877 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp" "-11.60386802 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp" "1.46449188 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "flux" "9.88474603 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy" "-21.41726768 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic" "21.53688965 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential" "-42.95415737 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy" "-5.02230991 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp" "-13.10442235 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp" "0.72390852 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "flux" "-12.59165665 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy" "-21.58100325 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic" "21.97625898 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential" "-43.55726217 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy" "-5.25263408 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp" "-14.40418877 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp" "1.65089413 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "flux" "-10.97538990 0.001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy" "-21.58034154 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic" "19.06945306 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential" "-40.64979462 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy" "-4.95910877 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp" "-11.60386802 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp" "1.46451581 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "flux" "9.88054883 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy" "-21.41383354 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic" "21.54062110 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential" "-42.95445468 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy" "-5.02235546 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp" "-13.10432010 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp" "0.72355452 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "flux" "-12.59721251 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy" "-21.57967929 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic" "21.97767070 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential" "-43.55734994 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy" "-5.25255007 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp" "-14.40431020 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp" "1.65084377 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "flux" "-10.97712288 0.001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy" "-21.65280108 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic" "19.03441381 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential" "-40.68721489 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy" "-5.28791536 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp" "-11.38769047 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp" "1.53971801 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "flux" "-10.19903930 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy" "-21.57285892 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic" "21.98319329 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential" "-43.55605221 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy" "-4.55367707 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp" "-14.68276787 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp" "1.23171979 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "flux" "-7.10804334 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy" "-21.52166382 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic" "19.23592476 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential" "-40.75758858 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy" "-5.38158923 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp" "-11.78892407 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp" "1.96425191 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "flux" "0.28132980 0.000001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy" "-21.65196333 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic" "19.03525227 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential" "-40.68721560 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy" "-5.28791536 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp" "-11.38769047 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp" "1.53971730 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "flux" "-10.20067791 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy" "-21.56976165 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic" "21.98663172 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential" "-43.55639337 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy" "-4.55372868 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp" "-14.68286398 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp" "1.23152635 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "flux" "-7.11570289 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy" "-21.52110513 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic" "19.23678572 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential" "-40.75789085 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy" "-5.38153225 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp" "-11.78962964 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp" "1.96459823 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "flux" "0.28095754 0.00001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_vmcbatch_dmcbatch_mwalkers.in.xml
    2
    2
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers.in.xml
    2
    2
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy" "-21.71827114 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic" "20.19340480 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential" "-41.91167602 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy" "-4.94686188 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp" "-13.34477165 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp" "1.93129114 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "flux" "-10.14767210 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy" "-21.98993631 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic" "21.53819515 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential" "-43.52813136 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy" "-5.77073014 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp" "-13.96799663 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp" "1.76192905 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "flux" "2.38005678 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy" "-21.98443887 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic" "21.24965809 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential" "-43.23409691 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy" "-5.54241465 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp" "-14.52283021 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp" "2.38248141 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "flux" "0.18637741 0.001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy" "-21.71738267 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic" "20.19428493 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential" "-41.91166755 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy" "-4.94686188 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp" "-13.34477165 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp" "1.93129961 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "flux" "-10.14942551 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy" "-21.99092404 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic" "21.53710746 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential" "-43.52803147 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy" "-5.77057594 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp" "-13.96787348 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp" "1.76175159 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "flux" "2.38866285 0.001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy" "-21.98218308 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic" "21.25165864 0.0005")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential" "-43.23384181 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy" "-5.54248445 0.00003")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion" "-25.55133363 0.000002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp" "-14.52254631 0.0002")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp" "2.38252241 0.0001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "flux" "0.17778549 0.001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy" "-21.68083037 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic" "20.16256817 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential" "-41.84339855 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy" "-4.90964385 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp" "-13.25550100 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp" "1.87307337 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "flux" "-8.99093561 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy" "-21.81825958 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic" "19.42117537 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential" "-41.23943495 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy" "-5.43148508 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp" "-11.91172090 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp" "1.65509810 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "flux" "1.80616014 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy" "-21.89208131 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic" "18.94622832 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential" "-40.83830963 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy" "-5.29808044 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp" "-12.58588514 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp" "2.59698314 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "flux" "-0.63233432 0.000001")
    else()
      # VMC without drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy" "-21.67973611 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic" "20.16365667 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential" "-41.84339279 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy" "-4.90964385 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp" "-13.25550100 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp" "1.87307913 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "flux" "-8.99299657 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples" "144.00000000 0.0")
      # VMC with drift
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy" "-21.81914511 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic" "19.42053111 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential" "-41.23967622 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy" "-5.43133155 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion" "-25.55132706 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp" "-11.91199362 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp" "1.65497602 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "flux" "1.81290367 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples" "144.00000000 0.0")
      # DMC
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy" "-21.89336143 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic" "18.94521744 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential" "-40.83857887 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy" "-5.29804220 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp" "-12.58642065 0.000001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp" "2.59721117 0.00001")
      list(APPEND DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "flux" "-0.63113857 0.00001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_vmcbatch_dmcbatch_mwalkers.in.xml
    2
    3
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_inverter-host_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_vmcbatch_dmcbatch_mwalkers_inverter_host.in.xml
    2
    3
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdbatch_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers.in.xml
    2
    3
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS # DMC
  )

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch-mwalkers_sdbatch_inverter-host_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers
    det_qmc_short_sdbatch_vmcbatch_dmcbatch_mwalkers_inverter_host.in.xml
    2
    3
    TRUE
    0
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS # VMC without drift
    1
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS # VMC with drift
    2
    DET_DIAMOND2_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "totenergy" "-20.82333418 0.00002128")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "kinetic" "16.28201421 0.00002563")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "potential" "-37.10534849 0.00000861")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "eeenergy" "-6.62968020 0.00000192")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "ionion" "-25.55133345 0.0000014")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "localecp" "-5.47510064 0.00000777")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "nonlocalecp" "0.55076580 0.00001229")
    else()
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "totenergy" "-20.82178481 0.00001859")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "kinetic" "16.28358306 0.00002343")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "potential" "-37.10536785 0.00001142")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "eeenergy" "-6.62966523 0.00000163")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "ionion" "-25.55133345 0.00000135")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "localecp" "-5.47512422 0.00001074")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "nonlocalecp" "0.55075504 0.00001354")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "totenergy" "-20.82333435 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "kinetic" "16.28200884 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "potential" "-37.10534319 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "eeenergy" "-6.62968027 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "localecp" "-5.47510116 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "nonlocalecp" "0.55076543 0.000001")
    else()
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "totenergy" "-20.82177682 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "kinetic" "16.28358399 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "potential" "-37.10536081 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "eeenergy" "-6.62966513 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "ionion" "-25.55132719 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "localecp" "-5.47512288 0.000001")
      list(APPEND DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS "nonlocalecp" "0.55075440 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-dmcbatch_sdj_excited
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmcbatch_dmcbatch_excited
    det_qmc_short_vmcbatch_dmcbatch_excited.in.xml
    1
    1
    TRUE
    1
    DET_DIAMOND2_DMC_BATCH_EXCITED_SCALARS # DMC
  )


  # CUDA references
else()
  # VMC runs
  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-21.89020466 0.000006")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.34475369 0.000004")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.23495832 0.000002")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-5.807169758 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-12.80697366 0.0000025")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.930512045 0.000002")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-10.84154438 0.000028")
    else()
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-21.88556609 0.00004857")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.34914975 0.000006")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.23473245 0.000009")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-5.807173593 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-12.80669499 0.000007")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.93046308 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-10.85016132 0.000001")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-21.89021689 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.34474147 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.23495836 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-5.80717002 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-12.80697507 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.93051366 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-10.84152599 0.000001")
    else()
      list(APPEND DET_DIAMOND2_SCALARS "totenergy" "-21.88560221 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "kinetic" "19.34913447 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "potential" "-41.23473668 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "eeenergy" "-5.80717246 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "localecp" "-12.80670171 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "nonlocalecp" "2.93046442 0.000001")
      list(APPEND DET_DIAMOND2_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_SCALARS "flux" "-10.85014255 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short
    det_qmc_short.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_SCALARS # VMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-19.35642193 0.000010")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "18.11199241 0.000011")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-37.46841434 0.0000012")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-5.239454858 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-7.801557933 0.000005")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.123925405 0.000006")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-7.266823563 0.000017")
    else()
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-19.33598072 0.000009")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "18.13250214 0.000015")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-37.46848286 0.000024")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-5.2395315 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-7.80151826 0.000016")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.123893845 0.000008")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-7.307580128 0.000039")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-19.35640958 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "18.11200421 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-37.46841379 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-5.23945567 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-7.80156231 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.12393113 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-7.26683500 0.000002")
    else()
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "totenergy" "-19.33596948 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "kinetic" "18.13251898 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "potential" "-37.46848846 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "eeenergy" "-5.23953030 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "localecp" "-7.80152053 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "nonlocalecp" "1.12388931 0.000001")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_EXCITED_SCALARS "flux" "-7.30760928 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-vmc_sdj_excited
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_excited
    det_qmc_short_excited.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_EXCITED_SCALARS # VMC
  )

  # VMC DELAYED-UPDATE
  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-21.7427256 0.00000365")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "17.06504779 0.00002395")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-38.80777339 0.00002067")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-7.303211623 0.000012")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-6.872058488 0.0000239")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "0.918823673 0.00000234")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-10.69690014 0.00009929")
    else()
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-21.74329658 0.000016")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "17.06563918 0.0000046")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-38.80893576 0.00001943")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-7.302218153 0.0000127")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-6.874538168 0.000038")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "0.919147508 0.00000571")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-10.69114954 0.0000894")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-21.74275370 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "17.06504267 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-38.80779637 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-7.30320156 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-6.87209644 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "0.91882857 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-10.69683005 0.000001")
    else()
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "totenergy" "-21.74333602 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "kinetic" "17.06562079 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "potential" "-38.80895681 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "eeenergy" "-7.30220276 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "localecp" "-6.87457977 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "nonlocalecp" "0.91915265 0.000001")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "samples" "9.00000000 0.0")
      list(APPEND DET_DIAMOND2_DELAY_SCALARS "flux" "-10.69100800 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-delayed_update-vmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short
    det_qmc_short_delayedupdate.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND2_DELAY_SCALARS # VMC
  )

  # DMC RUN
  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-19.40279161 0.00000399")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "18.56465394 0.0000031")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-37.96744555 0.000002")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-5.716501203 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-7.865690665 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "1.166073268 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-19.40331207 0.00000758")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "18.56412828 0.000008")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-37.96744035 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-5.71649912 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-7.865695195 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "1.16608091 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-19.40280208 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "18.56465081 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-37.96745289 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-5.71650065 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-7.86569882 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "1.16607351 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_SCALARS "totenergy" "-19.40332849 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "kinetic" "18.56411412 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "potential" "-37.96744262 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "eeenergy" "-5.71649975 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "localecp" "-7.86569564 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "nonlocalecp" "1.16607971 0.000001")
      list(APPEND DET_DIAMOND2_DMC_SCALARS "samples" "9.00000000 0.0")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-dmc_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmc_dmc
    det_qmc_short_vmc_dmc.in.xml
    1
    1
    TRUE
    1
    DET_DIAMOND2_DMC_SCALARS # DMC
  )

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-20.66429554 0.0000104")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "22.80989146 0.00000802")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.474187 0.0000032")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-0.24794429 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.58651551 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "-1.08840026 0.000003")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-20.65383644 0.00000879")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "22.82030478 0.000011")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.47414122 0.00000257")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-0.247937575 0.00000221")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55132695 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.58648575 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "-1.088390948 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    endif()
  else()
    if(QMC_COMPLEX)
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-20.66429132 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "22.80991044 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.47420176 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-0.24794375 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.58651611 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "-1.08841496 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    else()
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "totenergy" "-20.65380482 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "kinetic" "22.82032553 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "potential" "-43.47413035 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "eeenergy" "-0.24793700 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "ionion" "-25.55132694 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "localecp" "-16.58647423 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "nonlocalecp" "-1.08839219 0.000001")
      list(APPEND DET_DIAMOND2_DMC_EXCITED_SCALARS "samples" "9.00000000 0.0")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1_pp-dmc_sdj_excited
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
    det_qmc_short_vmc_dmc_excited
    det_qmc_short_vmc_dmc_excited.in.xml
    1
    1
    TRUE
    1
    DET_DIAMOND2_DMC_EXCITED_SCALARS # DMC
  )
endif()
