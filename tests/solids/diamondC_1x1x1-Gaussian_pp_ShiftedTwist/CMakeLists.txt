
# Reference VMC run in qmc_ref_long.s000.scalar.dat
# VMC short run with No Jastrow. Compare directly to Hartree Fock Energy of: E=-11.233944514046865
#Data generated using utils/make_ref_data.sh DIAMOND qmc_short.s000.scalar.dat ../qmc-ref/qmc_ref_long.s000.scalar.dat
IF(ENABLE_SOA)

LIST(APPEND LONG_DIAMOND_SCALARS "HartreeFock" "-11.233944514046865 0.0") 
LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-11.233944514046865 0.0")
#LIST(APPEND LONG_DIAMOND_SCALARS "variance" "4.489893400599 0.160802504189")
#LIST(APPEND LONG_DIAMOND_SCALARS "kinetic" "19.913500891733 0.039462938665")
#LIST(APPEND LONG_DIAMOND_SCALARS "potential" "-41.121249050943 0.040437159289")
#LIST(APPEND LONG_DIAMOND_SCALARS "eeenergy" "-4.465700045066 0.018685263294")
#LIST(APPEND LONG_DIAMOND_SCALARS "localecp" "-13.135076417429 0.057750691393")
#LIST(APPEND LONG_DIAMOND_SCALARS "nonlocalecp" "2.030854168240 0.034268032958")
#LIST(APPEND LONG_DIAMOND_SCALARS "ionion" "-25.551326757000 0.001000000000")
  QMC_RUN_AND_CHECK(long_diamondC_1x1x1-Gaussian_pp_ShiftedTwist
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_ShiftedTwist"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )
#
MESSAGE("Need To fix Phase!!!! Energy to find is -11.233944514046865 Ha" )
ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
