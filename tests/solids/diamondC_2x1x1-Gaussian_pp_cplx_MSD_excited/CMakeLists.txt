if(QMC_COMPLEX)

  #N_det             =          101
  #N_states          =            2
  #N_sop             =           17

  # * State            1
  # < S^2 >         =   1.432079152684324E-006
  # E               =   -21.3657227937065
  # Variance        =    1.00045224961504       +/-   0.000000000000000E+000
  # PT norm         =   0.461839003791580       +/-   0.000000000000000E+000
  # PT2             =  -0.404385553360503       +/-   0.000000000000000E+000
  # rPT2            =  -0.333295253748026       +/-   0.000000000000000E+000
  # E+PT2           =   -21.7671083470670       +/-   0.000000000000000E+000
  # E+rPT2          =   -21.6960180474545       +/-   0.000000000000000E+000
  #
  #qmc_ref_ground_long_noj  series 0
  #  LocalEnergy           =         -21.36523 +/-          0.00013
  #  Variance              =            4.1743 +/-           0.0045
  #  Kinetic               =          19.75794 +/-          0.00042
  #  LocalPotential        =         -41.12317 +/-          0.00044
  #  ElecElec              =          -4.60299 +/-          0.00016
  #  LocalECP              =         -13.07070 +/-          0.00056
  #  NonLocalECP           =           2.10185 +/-          0.00019
  #  IonIon                =            -25.55 +/-             0.00
  #  LocalEnergy_sq        =          460.6475 +/-           0.0069
  #  BlockWeight           =        1310720.00 +/-             0.00
  #  BlockCPU              =          22.54826 +/-          0.00013
  #  AcceptRatio           =         0.9355798 +/-        0.0000023
  #  Efficiency            =             23.55 +/-             0.00
  #  TotalTime             =           6471.35 +/-             0.00
  #  TotalSamples          =         376176640 +/-                0

  # * State            2
  # < S^2 >         =   6.026072064414621E-005
  # E               =   -20.7021278060674
  # Variance        =    1.05958506927153       +/-   0.000000000000000E+000
  # PT norm         =    1.56509030954092       +/-   0.000000000000000E+000
  ## PT2             =  -0.488098362967976       +/-   0.000000000000000E+000
  # rPT2            =  -0.141497978456392       +/-   0.000000000000000E+000
  # E+PT2           =   -21.1912261690353       +/-   0.000000000000000E+000
  # E+rPT2          =   -20.8446257845238       +/-   0.000000000000000E+000

  #qmc_ref_excited_long  series 0
  #  LocalEnergy           =         -20.70253 +/-          0.00016
  #  Variance              =             4.154 +/-            0.012
  #  Kinetic               =          20.63096 +/-          0.00043
  #  LocalPotential        =         -41.33350 +/-          0.00044
  #  ElecElec              =          -4.67055 +/-          0.00018
  #  LocalECP              =         -13.32539 +/-          0.00057
  #  NonLocalECP           =           2.21377 +/-          0.00022
  #  IonIon                =            -25.55 +/-             0.00
  #  LocalEnergy_sq        =           432.748 +/-            0.013
  #  BlockWeight           =        1310720.00 +/-             0.00
  #  BlockCPU              =          22.36911 +/-          0.00013
  #  AcceptRatio           =         0.9336608 +/-        0.0000020
  #  Efficiency            =             16.86 +/-             0.00
  #  TotalTime             =           6487.04 +/-             0.00
  #  TotalSamples          =         380108800 +/-                0

  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT vmc_short_ground.s000.scalar.dat  qmc_ref_ground_long_noj.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "totenergy" "-21.3657227937065 0.068654689630")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "variance" "4.173318575244 0.206016680533")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "kinetic" "19.757886856418 0.297513259260")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "potential" "-41.123130287178 0.282895938880")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "eeenergy" "-4.603005269151 0.085411994867")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "localecp" "-13.070694487044 0.441491994197")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "nonlocalecp" "2.101896531043 0.146195537980")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "ionion" "-25.551327062000 0.000100000000")

  qmc_run_and_check(
    short-diamondC_2x1x1-Gaussian_pp_cplx_ground_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"
    vmc_short_ground
    vmc_short_ground_noj.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS #VMC
  )

  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG vmc_long_ground_noj.s000.scalar.dat  qmc_ref_ground_long_noj.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "totenergy" "-21.3657227937065 0.025316995774")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "variance" "4.173318575244 0.129185833591")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "kinetic" "19.757886856418 0.073759934621")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "potential" "-41.123130287178 0.076620750856")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "eeenergy" "-4.603005269151 0.030347826563")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "localecp" "-13.070694487044 0.093758470517")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "nonlocalecp" "2.101896531043 0.033542969378")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "ionion" "-25.551327062000 0.000100000000")

  qmc_run_and_check(
    long-diamondC_2x1x1-Gaussian_pp_cplx_ground_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"
    vmc_long_ground_noj
    vmc_long_ground_noj.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS #VMC
  )

  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG vmc_long_excited_noj.s000.scalar.dat  qmc_ref_excited_long.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "totenergy" "-20.7021278060674 0.015464176265")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "variance" "4.153573777819 0.120610263740")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "kinetic" "20.630997778670 0.061617920662")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "potential" "-41.333533646609 0.061539681991")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "eeenergy" "-4.670554461220 0.018805728911")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "localecp" "-13.325392368394 0.084934045045")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "nonlocalecp" "2.213740245050 0.021921851285")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "ionion" "-25.551327062000 0.000100000000")

  qmc_run_and_check(
    long-diamondC_2x1x1-Gaussian_pp_cplx_excited_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"
    vmc_long_excited_noj
    vmc_long_excited_noj.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS #VMC
  )

  #Data generated using utils/make_ref_data.sh diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT vmc_short_excited_noj.s000.scalar.dat  qmc_ref_excited_long.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "totenergy" "-20.7021278060674 0.053373569209")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "variance" "4.153573777819 0.215901996866")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "kinetic" "20.630997778670 0.186121488914")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "potential" "-41.333533646609 0.197753128678")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "eeenergy" "-4.670554461220 0.057273234412")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "localecp" "-13.325392368394 0.264864425798")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "nonlocalecp" "2.213740245050 0.089172131789")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "ionion" "-25.551327062000 0.0000100000000")

  qmc_run_and_check(
    short-diamondC_2x1x1-Gaussian_pp_cplx_excited_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"
    vmc_short_excited_noj
    vmc_short_excited_noj.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS #VMC
  )

else()
  message(VERBOSE "Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
endif()
