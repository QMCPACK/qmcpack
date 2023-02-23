if(QMC_COMPLEX)
  if(NOT ENABLE_OFFLOAD)
    list(APPEND BCC_MO_VMC_SCALARS "totenergy" "-136.52078308 0.016914")
    list(APPEND BCC_MO_VMC_SCALARS "kinetic" "47.99405288 0.094614")
    list(APPEND BCC_MO_VMC_SCALARS "potential" "-184.51483596 0.095482")
    list(APPEND BCC_MO_VMC_SCALARS "eeenergy" "1.57730687 0.047504")
    list(APPEND BCC_MO_VMC_SCALARS "ionion" "-119.82754015 0.00002")
    list(APPEND BCC_MO_VMC_SCALARS "localecp" "-91.69929753 0.183567")
    list(APPEND BCC_MO_VMC_SCALARS "nonlocalecp" "25.45993896 0.126237")
    list(APPEND BCC_MO_VMC_SCALARS "soecp" "-0.02524410 0.000991")
    list(APPEND BCC_MO_VMC_SCALARS "samples" "40000 0.0")

    qmc_run_and_check(
      short-bccMo_2x1x1_SOREP-vmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/bccMo_2x1x1_SOREP"
      Mo-vmc-short
      Mo-vmc-short.in.xml
      1
      16
      TRUE
      0
      BCC_MO_VMC_SCALARS # VMC
    )

    #
    # Long tests
    #

    list(APPEND LONG_BCC_MO_VMC_SCALARS "totenergy" "-136.52078308 0.002674")
    list(APPEND LONG_BCC_MO_VMC_SCALARS "soecp" "-0.02524410 0.000157")
    list(APPEND LONG_BCC_MO_VMC_SCALARS "samples" "1600000 0.0")
    qmc_run_and_check(
      long-bccMo_2x1x1_SOREP-vmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/bccMo_2x1x1_SOREP"
      Mo-vmc-long
      Mo-vmc-long.in.xml
      1
      16
      TRUE
      0
      LONG_BCC_MO_VMC_SCALARS # VMC
    )

    list(APPEND BCC_MO_DMC_SCALARS "totenergy" "-136.77956168 0.063003")
    list(APPEND BCC_MO_DMC_SCALARS "kinetic" "48.07180161 0.451909")
    list(APPEND BCC_MO_DMC_SCALARS "potential" "-184.85136330 0.454768")
    list(APPEND BCC_MO_DMC_SCALARS "eeenergy" "1.53598049 0.291135")
    list(APPEND BCC_MO_DMC_SCALARS "ionion" "-119.82754015 0.00002")
    list(APPEND BCC_MO_DMC_SCALARS "localecp" "-92.14260260 1.147154")
    list(APPEND BCC_MO_DMC_SCALARS "nonlocalecp" "25.60831010 0.643639")
    list(APPEND BCC_MO_DMC_SCALARS "soecp" "-0.02551113 0.003011")

    qmc_run_and_check(
      short-bccMo_2x1x1_SOREP-dmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/bccMo_2x1x1_SOREP"
      Mo-vmc-dmc-short
      Mo-vmc-dmc-short.in.xml
      1
      16
      TRUE
      1
      BCC_MO_DMC_SCALARS # DMC
    )

    list(APPEND LONG_BCC_MO_DMC_SCALARS "totenergy" "-136.77956168 0.009961")
    list(APPEND LONG_BCC_MO_DMC_SCALARS "soecp" "-0.02551113 0.000476")
    qmc_run_and_check(
      long-bccMo_2x1x1_SOREP-dmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/bccMo_2x1x1_SOREP"
      Mo-vmc-dmc-long
      Mo-vmc-dmc-long.in.xml
      1
      16
      TRUE
      1
      LONG_BCC_MO_DMC_SCALARS # DMC
    )

    #Deterministic test

    if(NOT QMC_MIXED_PRECISION)
      list(APPEND DET_BCC_MO_VMC_SCALARS "totenergy" "-136.88580235 0.000002")
      list(APPEND DET_BCC_MO_VMC_SCALARS "kinetic" "48.19451127 0.000001")
      list(APPEND DET_BCC_MO_VMC_SCALARS "potential" "-185.08031362 0.000002")
      list(APPEND DET_BCC_MO_VMC_SCALARS "eeenergy" "0.73302788 0.000001")
      list(APPEND DET_BCC_MO_VMC_SCALARS "ionion" "-119.82754227 0.000002")
      list(APPEND DET_BCC_MO_VMC_SCALARS "localecp" "-90.13681859 0.000001")
      list(APPEND DET_BCC_MO_VMC_SCALARS "nonlocalecp" "24.16783601 0.000001")
      list(APPEND DET_BCC_MO_VMC_SCALARS "soecp" "-0.01681613 0.000001")
      list(APPEND DET_BCC_MO_VMC_SCALARS "samples" "144 0.0")

      qmc_run_and_check(
        deterministic-bccMo_2x1x1_SOREP-vmc_sdj
        "${qmcpack_SOURCE_DIR}/tests/solids/bccMo_2x1x1_SOREP"
        det_Mo-vmc
        det_Mo-vmc.in.xml
        1
        1
        TRUE
        0
        DET_BCC_MO_VMC_SCALARS # VMC
      )

      list(APPEND DET_BCC_MO_DMC_SCALARS "totenergy" "-135.34032696 0.000002")
      list(APPEND DET_BCC_MO_DMC_SCALARS "kinetic" "29.19088964 0.000001")
      list(APPEND DET_BCC_MO_DMC_SCALARS "potential" "-164.53121661 0.000002")
      list(APPEND DET_BCC_MO_DMC_SCALARS "eeenergy" "-2.19402623 0.000001")
      list(APPEND DET_BCC_MO_DMC_SCALARS "ionion" "-119.82754227 0.000002")
      list(APPEND DET_BCC_MO_DMC_SCALARS "localecp" "-64.57264081 0.000001")
      list(APPEND DET_BCC_MO_DMC_SCALARS "nonlocalecp" "22.07409626 0.000001")
      list(APPEND DET_BCC_MO_DMC_SCALARS "soecp" "-0.01110304 0.000001")

      qmc_run_and_check(
        deterministic-bccMo_2x1x1_SOREP-dmc_sdj
        "${qmcpack_SOURCE_DIR}/tests/solids/bccMo_2x1x1_SOREP"
        det_Mo-vmc-dmc
        det_Mo-vmc-dmc.in.xml
        1
        1
        TRUE
        1
        DET_BCC_MO_DMC_SCALARS # DMC
      )
    endif()
  endif()
endif()
