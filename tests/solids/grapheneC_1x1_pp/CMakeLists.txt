#  GCC 6.2 References
#  LIST(APPEND DET_GRAPHENE_SCALARS "ionion" "14.00303308 0.000001")

# Clang 7.0 References
#  LIST(APPEND DET_GRAPHENE_SCALARS "ionion" "14.00303155 0.000001")   

# Intel 18 References
#  LIST(APPEND DET_GRAPHENE_SCALARS "ionion" "14.00302796 0.000001")

#  Q.E. reference
  LIST(APPEND DET_GRAPHENE_SCALARS "ionion" "14.002118125 0.000001") 

  QMC_RUN_AND_CHECK(issue-deterministic-grapheneC_1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/grapheneC_1x1_pp"
                    det_qmc_short
                    det_qmc_short.in.xml
                    1 1
                    TRUE
                    0 DET_GRAPHENE_SCALARS # VMC
                    )

#  Intel 18 References
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "totenergy" "-11.14884707 0.000001")
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "kinetic" "10.93205382 0.000001")
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "potential" "-22.08090088 0.000001")
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "eeenergy" "51.43836209 0.000001")
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "ionion" "47.77174570 0.000001") # DFT reference : 47.80404842
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "localecp" "-122.58517139 0.000001")
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "nonlocalecp" "1.29416271 0.000001")
#  LIST(APPEND DET_GRAPHENE_VACUUM_SCALARS "samples" "9 0.0")

#  QMC_RUN_AND_CHECK(deterministic-grapheneC_1x1_pp-vmc_sdj_vacuum
#                    "${CMAKE_SOURCE_DIR}/tests/solids/grapheneC_1x1_pp"
#                    det_qmc_short_vacuum
#                    det_qmc_short_vacuum.in.xml
#                    1 1
#                    TRUE
#                    0 DET_GRAPHENE_VACUUM_SCALARS # VMC
#                    )
