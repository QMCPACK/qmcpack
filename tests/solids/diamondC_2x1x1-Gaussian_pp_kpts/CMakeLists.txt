
 IF (NOT QMC_CUDA)
   IF (QMC_COMPLEX)

# Reference VMC run in qmc_ref_long.s000.scalar.dat
# VMC short run with No Jastrow. Compare directly to Hartree Fock Energy of: E=-21.20673553792076
#Data generated using utils/make_ref_data.sh DIAMOND qmc_short.s000.scalar.dat ../qmc-ref/qmc_ref_long.s000.scalar.dat
LIST(APPEND DIAMOND_SCALARS "HartreeFock" "-21.20673553792076 0.0")
LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.207748159198 0.034802492057")
LIST(APPEND DIAMOND_SCALARS "variance" "4.489893400599 0.215879917849")
LIST(APPEND DIAMOND_SCALARS "kinetic" "19.913500891733 0.146083699450")
LIST(APPEND DIAMOND_SCALARS "potential" "-41.121249050943 0.162975074780")
LIST(APPEND DIAMOND_SCALARS "eeenergy" "-4.465700045066 0.036509787281")
LIST(APPEND DIAMOND_SCALARS "localecp" "-13.135076417429 0.178033617003")
LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "2.030854168240 0.050570571341")
LIST(APPEND DIAMOND_SCALARS "ionion" "-25.551326757000 0.001000000000")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp_kpts-vmc_gaussian_sdj
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

#
# Long tests
#
# VMC long run with No Jastrow. Compare directly to Hartree Fock Energy of: E=-21.20673553792076
#Data generated using utils/make_ref_data.sh DIAMOND qmc_long.s000.scalar.dat ../qmc-ref/qmc_ref_long.s000.scalar.dat


LIST(APPEND LONG_DIAMOND_SCALARS "HartreeFock" "-21.20673553792076 0.0") 
LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.207748159198 0.017379145929")
LIST(APPEND LONG_DIAMOND_SCALARS "variance" "4.489893400599 0.160802504189")
LIST(APPEND LONG_DIAMOND_SCALARS "kinetic" "19.913500891733 0.039462938665")
LIST(APPEND LONG_DIAMOND_SCALARS "potential" "-41.121249050943 0.040437159289")
LIST(APPEND LONG_DIAMOND_SCALARS "eeenergy" "-4.465700045066 0.018685263294")
LIST(APPEND LONG_DIAMOND_SCALARS "localecp" "-13.135076417429 0.057750691393")
LIST(APPEND LONG_DIAMOND_SCALARS "nonlocalecp" "2.030854168240 0.034268032958")
LIST(APPEND LONG_DIAMOND_SCALARS "ionion" "-25.551326757000 0.001000000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-vmc_gaussian_sdj
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )
#
#  DMC Result with 1, 2 and 3Body Jastrow
#
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.834535002421 0.011309444536")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "variance" "1.304577591187 0.017886443257")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "kinetic" "20.405490229411 0.085631287654")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "potential" "-42.240025231827 0.086355601969")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "eeenergy" "-5.450656130373 0.023207167719")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "localecp" "-13.228959467675 0.114648569345")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "nonlocalecp" "1.990917122865 0.045023816102")
LIST(APPEND LONG_DIAMOND_DMC_SCALARS "ionion" "-25.551326757000 0.001000000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-dmc_gaussian_sdj
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )



    ELSE()
      MESSAGE_VERBOSE("Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
    ENDIF()
  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
