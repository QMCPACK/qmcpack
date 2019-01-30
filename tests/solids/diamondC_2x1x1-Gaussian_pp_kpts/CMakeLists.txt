
# Reference VMC run in qmc_ref_long.s000.scalar.dat
IF(ENABLE_SOA)
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.710528403484 0.014941420854")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "20.482953080839 0.115723306504")
  LIST(APPEND DIAMOND_SCALARS "potential" "-42.193481484258 0.123056258811")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-5.356726989742 0.031142760274")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-13.125425643581 0.116081264579")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "1.839997905661 0.042100003374")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-25.551326757000 0.000100000000")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp_kpts-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

#
# Long tests
#

# Reference VMC run in qmc_ref_long_vmc_dmc.s001.scalar.dat"
LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.710207226000 0.006931462353")
LIST(APPEND LONG_DIAMOND_SCALARS "kinetic" "20.499369979000 0.063213462871")
LIST(APPEND LONG_DIAMOND_SCALARS "potential" "-42.209577205000 0.059854365950")
LIST(APPEND LONG_DIAMOND_SCALARS "eeenergy" "-5.356961351300 0.020388384354")
LIST(APPEND LONG_DIAMOND_SCALARS "localecp" "-13.132922225000 0.079866755274")
LIST(APPEND LONG_DIAMOND_SCALARS "nonlocalecp" "1.831633128100 0.018918103422")
LIST(APPEND LONG_DIAMOND_SCALARS "ionion" "-25.551326757000 0.000100000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )
# generated scalars by
#../../../../../utils/make_ref_data.sh LONG_DIAMOND_DMC qmc_long_vmc_dmc.s001.scalar.dat qmc-ref/qmc_ref_long_vmc_dmc.s001.scalar.dat
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "kinetic" "20.538888744528 0.083505174054")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "potential" "-42.383282722862 0.082006395952")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "eeenergy" "-5.452768166889 0.049835970216")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "localecp" "-13.212171620553 0.079056548891")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "nonlocalecp" "1.832983821228 0.028452375262")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "ionion" "-25.551326757000 0.000100000000")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.84416881 0.00083760")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
