
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

# ../../../../../utils/make_ref_data.sh LONG_DIAMOND qmc_long.s000.scalar.dat qmc-ref/qmc_ref_long.s000.scalar.dat
LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.710283757096 0.007003082789")
LIST(APPEND LONG_DIAMOND_SCALARS "kinetic" "20.484766989575 0.040286069522")
LIST(APPEND LONG_DIAMOND_SCALARS "potential" "-42.195050746630 0.040257715256")
LIST(APPEND LONG_DIAMOND_SCALARS "eeenergy" "-5.355965775116 0.011238545502")
LIST(APPEND LONG_DIAMOND_SCALARS "localecp" "-13.128878317685 0.043236660046")
LIST(APPEND LONG_DIAMOND_SCALARS "nonlocalecp" "1.841120102814 0.021549673940")
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

LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.844391706214 0.007790489092")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "kinetic" "20.538585298968 0.060699060976")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "potential" "-42.382977005143 0.047872875660")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "eeenergy" "-5.452921721419 0.017452521845")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "localecp" "-13.211238602636 0.054588859799")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "nonlocalecp" "1.832510075558 0.038600443160")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "ionion" "-25.551326757000 0.000100000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
