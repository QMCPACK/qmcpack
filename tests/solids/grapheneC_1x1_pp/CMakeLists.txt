# Q.E. reference
#  LIST(APPEND DET_GRAPHENE_Z10_SCALARS "ionion" "14.002118125 0.000001") 

  LIST(APPEND DET_GRAPHENE_Z10_SCALARS "ionion" "14.00189436 0.000001")

  QMC_RUN_AND_CHECK(deterministic-grapheneC_1x1_pp-vmc_sdj_z10
                    "${CMAKE_SOURCE_DIR}/tests/solids/grapheneC_1x1_pp"
                    det_qmc_short_z10
                    det_qmc_short_z10.in.xml
                    1 1
                    TRUE
                    0 DET_GRAPHENE_Z10_SCALARS # VMC
                    )

# Q.E. References
#  LIST(APPEND DET_GRAPHENE_Z30_SCALARS "ionion" "81.60597872 0.000001")

  LIST(APPEND DET_GRAPHENE_Z30_SCALARS "ionion" "81.60576396 0.000001")

  QMC_RUN_AND_CHECK(deterministic-grapheneC_1x1_pp-vmc_sdj_z30
                    "${CMAKE_SOURCE_DIR}/tests/solids/grapheneC_1x1_pp"
                    det_qmc_short_z30
                    det_qmc_short_z30.in.xml
                    1 1
                    TRUE
                    0 DET_GRAPHENE_Z30_SCALARS # VMC
                    )
