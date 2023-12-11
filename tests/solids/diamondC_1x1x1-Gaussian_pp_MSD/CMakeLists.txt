# VMC long run with no Jastrows
# Compare directly to sCI Energy of: E=-10.39616312, T=11.2105174, and V=-21.60668054
# The other references are from VMC reference runs.
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "totenergy" "-10.39616312 0.0080035848")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "kinetic" "11.2105174 0.0344527216")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "potential" "-21.60668054 0.034737248")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "eeenergy" "-2.51687510 0.0123828708")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "localecp" "-7.13621881 0.0461973356")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "nonlocalecp" "0.82585828 0.015849296")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "ionion" "-12.77566743 0.000001")
qmc_run_and_check(
  long-diamondC_1x1x1-Gaussian_pp_MSD
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_MSD"
  vmc_long
  vmc_long.in.xml
  1
  16
  TRUE
  0
  diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS #VMC
)

# VMC short run with no Jastrows
# Compare directly to sCI Energy of: E=-10.39616312, T=11.2105174, and V=-21.60668054
# The other references are from VMC reference runs.
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "totenergy" "-10.39616312 0.025291327968")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "kinetic" "11.2105174 0.108870600256")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "potential" "-21.60668054 0.10976970368")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "eeenergy" "-2.51687510 0.039129871728")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "localecp" "-7.13621881 0.145983580496")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "nonlocalecp" "0.82585828 0.05008377536")
list(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "ionion" "-12.77566743 0.000001")
qmc_run_and_check(
  short-diamondC_1x1x1-Gaussian_pp_MSD
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_MSD"
  vmc_short
  vmc_short.in.xml
  1
  16
  TRUE
  0
  diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS #VMC
)

#Deterministic test
if(QMC_MIXED_PRECISION)
  #VMC
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-10.62293537 0.00002")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "12.00640791 0.00002")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-22.62934328 0.00001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-2.50996032 0.000002")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-7.67975015 0.00001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "0.33604476 0.000005")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-12.77567016 0.000001")
  if(QMC_COMPLEX)
    #DMC
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.99915119 0.00002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "11.94458085 0.00002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.94378393 0.000003")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.97216487 0.000002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-6.84338204 0.000002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "-0.35256283 0.000002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77567016 0.000001")
  else()
    #DMC
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.63742696 0.00002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "12.07373931 0.00002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.71116627 0.00004")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.72638743 0.000002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-7.91852665 0.00002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "0.70941796 0.00002")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77567016 0.000001")
  endif()
else()
  #VMC
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-10.62292565 0.000001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "12.00640819 0.000001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-22.62933384 0.000001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-2.50995954 0.000001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-7.67975160 0.000001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "0.33604476 0.000001")
  list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-12.77566743 0.000001")
  if(QMC_COMPLEX)
    #DMC
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.99920259 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "11.94457971 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.94378230 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.97216524 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-6.84338340 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "-0.35256619 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77566743 0.000001")
  else()
    #DMC
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.63742902 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "12.07372470 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.71115372 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.72638698 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-7.91851937 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "0.70942010 0.000001")
    list(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77566743 0.000001")
  endif()
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1-Gaussian_pp_MSD
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_MSD"
  det_vmc_dmc
  det_vmc_dmc.in.xml
  1
  1
  TRUE
  0
  det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS #VMC
  1
  det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS #DMC
)
