list(APPEND DIAMOND_SCALARS "totenergy" "-10.491445 0.003")
list(APPEND DIAMOND_SCALARS "kinetic" "11.434392 0.018")
list(APPEND DIAMOND_SCALARS "potential" "-21.925822 0.018")
list(APPEND DIAMOND_SCALARS "eeenergy" "-2.700534 0.0043")
list(APPEND DIAMOND_SCALARS "ionion" "-12.77567000 0.001")
list(APPEND DIAMOND_SCALARS "localecp" "-7.046740 0.030")
list(APPEND DIAMOND_SCALARS "nonlocalecp" "0.597119 0.0056")
list(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
list(APPEND DIAMOND_SCALARS "mpc" "-2.453044 0.004431")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

qmc_run_and_check(
  short-diamondC_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC
)

# Hybridrep is not implemented in legacy CUDA but should be correctly error trapped
qmc_run_and_check(
  short-diamondC_1x1x1_hybridrep_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short
  qmc_short_hybridrep.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC
)

qmc_run_and_check(
  short-diamondC_1x1x1_pp-vmc_sdj-meshf
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short-meshf
  qmc_short-meshf.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC
)

if(NOT ENABLE_OFFLOAD)
  list(APPEND DIAMOND_KSPACE_SCALARS "totenergy" "-10.500719 0.001769")
  list(APPEND DIAMOND_KSPACE_SCALARS "variance" "0.312264 0.028662")
  qmc_run_and_check(
    short-diamondC_1x1x1_pp-vmc_sdj_kspace
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_short_kspace
    qmc_short_kspace.in.xml
    1
    16
    TRUE
    0
    DIAMOND_KSPACE_SCALARS # VMC
  )

  qmc_run_and_check(
    short-diamondC_1x1x1_pp-vmc_sdj_kspace
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_short_kspace_4_4
    qmc_short_kspace_4_4.in.xml
    4
    4
    TRUE
    0
    DIAMOND_KSPACE_SCALARS # VMC
  )
else()
  message(VERBOSE "Skipping k-space Jastrow tests because they are not supported by OFFLOAD build (ENABLE_OFFLOAD=1)")
endif()

# Reference OPT run in qmc-ref
list(APPEND DIAMOND_OPT_SCALARS "totenergy" "-10.49370 0.0024")
list(APPEND DIAMOND_OPT_SCALARS "kinetic" "11.5504 0.021")
list(APPEND DIAMOND_OPT_SCALARS "localecp" "-7.1979 0.023")
list(APPEND DIAMOND_OPT_SCALARS "nonlocalecp" "0.61540 0.0068")
list(APPEND DIAMOND_OPT_SCALARS "samples" "128000 0.0")

qmc_run_and_check(
  short-diamondC_1x1x1_pp-opt_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short
  qmc_short_opt.in.xml
  1
  16
  ${SUCCESS_STATUS_MP}
  3
  DIAMOND_OPT_SCALARS # VMC
)

qmc_run_and_check(
  short-diamondC_1x1x1_pp-optbatch_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short
  qmc_short_optbatch.in.xml
  1
  3
  ${SUCCESS_STATUS_MP}
  3
  DIAMOND_OPT_SCALARS # VMC
)

# Reference DMC run in qmc-ref "-10.531583 0.000265"
list(APPEND DIAMOND_DMC_SCALARS "totenergy" "-10.5316 0.0024")
list(APPEND DIAMOND_DMC_SCALARS "kinetic" "11.4857 0.0231")
list(APPEND DIAMOND_DMC_SCALARS "potential" "-22.0170 0.0239")
list(APPEND DIAMOND_DMC_SCALARS "localecp" "-7.1518 0.0323")
list(APPEND DIAMOND_DMC_SCALARS "nonlocalecp" "0.62688 0.0080")
list(APPEND DIAMOND_DMC_SCALARS "eeenergy" "-2.71641 0.0073")
list(APPEND DIAMOND_DMC_SCALARS "mpc" "-2.47994 0.00779")

qmc_run_and_check(
  short-diamondC_1x1x1_pp-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short_vmc_dmc
  qmc_short_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  DIAMOND_DMC_SCALARS # DMC
)

# Excited state
list(APPEND DIAMOND_EXCITED_SCALARS "totenergy" "-10.28757413 0.001966")
list(APPEND DIAMOND_EXCITED_SCALARS "kinetic" "11.41769866 0.01674")
list(APPEND DIAMOND_EXCITED_SCALARS "potential" "-21.70527280 0.016756")
list(APPEND DIAMOND_EXCITED_SCALARS "eeenergy" "-2.78729624 0.004153")
list(APPEND DIAMOND_EXCITED_SCALARS "ionion" "-12.77566741 0.00001")
list(APPEND DIAMOND_EXCITED_SCALARS "localecp" "-6.74474055 0.018968")
list(APPEND DIAMOND_EXCITED_SCALARS "nonlocalecp" "0.60243140 0.005574")
list(APPEND DIAMOND_EXCITED_SCALARS "samples" "128000 0.0")
list(APPEND DIAMOND_EXCITED_SCALARS "mpc" "-2.52740814 0.004392")

qmc_run_and_check(
  short-diamondC_1x1x1_pp-vmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short_excited
  qmc_short_excited.in.xml
  1
  16
  TRUE
  0
  DIAMOND_EXCITED_SCALARS # VMC
)

list(APPEND DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-10.33047547 0.002606")
list(APPEND DIAMOND_DMC_EXCITED_SCALARS "kinetic" "11.48348626 0.026185")
list(APPEND DIAMOND_DMC_EXCITED_SCALARS "potential" "-21.81396173 0.02679")
list(APPEND DIAMOND_DMC_EXCITED_SCALARS "localecp" "-6.87151058 0.033201")
list(APPEND DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "0.63492715 0.008002")
list(APPEND DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-2.80171089 0.007489")
list(APPEND DIAMOND_DMC_EXCITED_SCALARS "mpc" "-2.55372811 0.00795")

qmc_run_and_check(
  short-diamondC_1x1x1_pp-dmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_short_vmc_dmc_excited
  qmc_short_vmc_dmc_excited.in.xml
  1
  16
  TRUE
  1
  DIAMOND_DMC_EXCITED_SCALARS # DMC
)

qmc_run_and_check(
  short-diamondC_1x1x1_pp-vmc-dmc-allp_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_vmc_dmc_allp
  qmc_short_vmc_dmc_allp.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC without drift
  1
  DIAMOND_SCALARS # VMC with drift
  2
  DIAMOND_DMC_SCALARS # DMC
)

# Coverage test - shorter version of dmc test
coverage_run("diamondC_1x1x1_pp-vmc_sdj-1-1" "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 1 1
             qmc_short_vmc_dmc_coverage.in.xml)

# Test time limit
cpu_limit_run("diamondC_1x1x1_pp-vmc-1-1" "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 1 1 180
              qmc_cpu_limit_vmc.xml)

cpu_limit_run("diamondC_1x1x1_pp-vmc-4-4" "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 4 4 120
              qmc_cpu_limit_vmc.xml)

cpu_limit_run("diamondC_1x1x1_pp-dmc-1-1" "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 1 1 240
              qmc_cpu_limit_dmc.xml)

cpu_limit_run("diamondC_1x1x1_pp-dmc-4-4" "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 4 4 120
              qmc_cpu_limit_dmc.xml)

#
# Long tests
#

list(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.491445 0.000065")
list(APPEND LONG_DIAMOND_SCALARS "samples" "122880000 0.0")
list(APPEND LONG_DIAMOND_SCALARS "flux" "0.0 0.03")
list(APPEND LONG_DIAMOND_SCALARS "mpc" "-2.453044 0.000143")

qmc_run_and_check(
  long-diamondC_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND_SCALARS # VMC
)

qmc_run_and_check(
  long-diamondC_1x1x1_pp-vmc_sdj-meshf
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_long-meshf
  qmc_long-meshf.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND_SCALARS # VMC
)

# VMC ref
list(APPEND LONG_DIAMOND_ALLP_SCALARS "totenergy" "-10.491445 0.0065")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "kinetic" "11.434392 0.006")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "potential" "-21.925822 0.006")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "eeenergy" "-2.700534 0.0015")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "ionion" "-12.77567000 0.0004")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "localecp" "-7.046740 0.0065")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "nonlocalecp" "0.597119 0.0019")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "mpc" "-2.453044 0.0014")
list(APPEND LONG_DIAMOND_ALLP_SCALARS "samples" "1280000 0.0")
#DMC ref
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "totenergy" "-10.5316 0.0008")
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "kinetic" "11.4857 0.0077")
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "potential" "-22.0170 0.0080")
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "localecp" "-7.1518 0.0108")
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "nonlocalecp" "0.62688 0.0027")
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "eeenergy" "-2.71641 0.0028")
list(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "mpc" "-2.47943 0.00268")

qmc_run_and_check(
  long-diamondC_1x1x1_pp-vmc-dmc-allp_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_vmc_dmc_allp
  qmc_long_vmc_dmc_allp.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND_ALLP_SCALARS # VMC without drift
  1
  LONG_DIAMOND_ALLP_SCALARS # VMC with drift
  2
  LONG_DIAMOND_DMC_ALLP_SCALARS # DMC
)

if(NOT ENABLE_OFFLOAD)
  list(APPEND LONG_DIAMOND_KSPACE_SCALARS "totenergy" "-10.500719 0.001769")
  list(APPEND LONG_DIAMOND_KSPACE_SCALARS "variance" "0.312264 0.028662")
  qmc_run_and_check(
    long-diamondC_1x1x1_pp-vmc_sdj_kspace
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_long_kspace
    qmc_long_kspace.in.xml
    1
    16
    TRUE
    0
    LONG_DIAMOND_KSPACE_SCALARS # VMC
  )
else()
  message(VERBOSE "Skipping k-space Jastrow tests because they are not supported by OFFLOAD build (ENABLE_OFFLOAD=1)")
endif()

# Reference DMC run in qmc-ref "-10.531583 0.000265"
list(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.531583 0.000815")
list(APPEND LONG_DIAMOND_DMC_SCALARS "mpc" "-2.47994 0.00246")

qmc_run_and_check(
  long-diamondC_1x1x1_pp-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_long_vmc_dmc
  qmc_long_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND_DMC_SCALARS # DMC
)

# Excited state
list(APPEND LONG_DIAMOND_EXCITED_SCALARS "totenergy" "-10.28757413 0.000067")
list(APPEND LONG_DIAMOND_EXCITED_SCALARS "samples" "122880000 0.0")
list(APPEND LONG_DIAMOND_EXCITED_SCALARS "flux" "-0.01538267 0.03786")

qmc_run_and_check(
  long-diamondC_1x1x1_pp-vmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_long_excited
  qmc_long_excited.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND_EXCITED_SCALARS # VMC
)

list(APPEND LONG_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-10.33047547 0.00086")

qmc_run_and_check(
  long-diamondC_1x1x1_pp-dmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  qmc_long_vmc_dmc_excited
  qmc_long_vmc_dmc_excited.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND_DMC_EXCITED_SCALARS # DMC
)

# estimator tests
include("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
check_python_reqs("numpy;h5py" diamond-estimator add_estimator_tests)
set(NMPI 4)
set(NOMP 4)
set(IFEXT "")

if(NOT QMC_COMPLEX)
  set(OFEXT "_real")
else()
  set(OFEXT "_comp")
endif()
if(add_estimator_tests)
  # density tests
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    density
    -e
    20
    -c
    8
    -p
    qmc_dens_short
    -r
    qmc-ref/qmc_dens_short.s000.stat_ref_density.dat)

  simple_run_and_check(
    short-diamondC_1x1x1_pp-dmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_dmc_short.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    density
    -e
    20
    -c
    8
    -p
    qmc_dens_dmc_short
    -r
    qmc-ref/qmc_dens_dmc_short.s001.stat_ref_density.dat)

  simple_run_and_check(
    short-diamondC_1x1x1_pp-rmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_rmc_short.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    density
    -e
    20
    -c
    8
    -p
    qmc_dens_rmc_short
    -r
    qmc-ref/qmc_dens_rmc_short.s001.stat_ref_density.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-vmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_long${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    density
    -e
    20
    -c
    8
    -p
    qmc_dens_long
    -r
    qmc-ref/qmc_dens_long.s000.stat_ref_density.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-dmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_dmc_long.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    density
    -e
    20
    -c
    8
    -p
    qmc_dens_dmc_long
    -r
    qmc-ref/qmc_dens_dmc_long.s001.stat_ref_density.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-rmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_rmc_long.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    density
    -e
    20
    -c
    8
    -p
    qmc_dens_rmc_long
    -r
    qmc-ref/qmc_dens_rmc_long.s001.stat_ref_density.dat)

  # spindensity tests
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    spindensity
    -e
    20
    -c
    8
    -p
    qmc_spindens_short
    -r
    qmc-ref/qmc_spindens_short.s000.stat_ref_spindensity.dat)

  # Still some issue with SpinDensityNew in the complete application
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmcbatch-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_vmcbatch_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    spindensity
    -e
    20
    -c
    8
    -p
    qmc_spindens_short
    -r
    qmc-ref/qmc_spindens_short.s000.stat_ref_spindensity.dat)

  # This is the new 1RDM test, it passes with 16 sigma, which is better than the
  # old test but we need to fix these in general.
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmcbatch-estimator-onebodydensitymatrices
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_onebodydensitymatrices_vmcbatch_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    obdm
    -e
    20
    -n
    16
    -c
    8
    -p
    qmc_onebodydensitymatrices_short
    -r
    qmc-ref/qmc_1rdm_noJ_short${OFEXT}.s000.stat_ref_1rdm.dat)

  simple_run_and_check(
    short-diamondC_1x1x1_pp-dmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_dmc_short.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    2
    -c
    8
    -p
    qmc_spindens_dmc_short
    -r
    qmc-ref/qmc_spindens_dmc_short.s001.stat_ref_spindensity.dat)

  simple_run_and_check(
    short-diamondC_1x1x1_pp-rmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_rmc_short.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    20
    -c
    8
    -p
    qmc_spindens_rmc_short
    -r
    qmc-ref/qmc_spindens_rmc_short.s001.stat_ref_spindensity.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-vmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_long${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    spindensity
    -e
    20
    -c
    8
    -p
    qmc_spindens_long
    -r
    qmc-ref/qmc_spindens_long.s000.stat_ref_spindensity.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-dmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_dmc_long.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    20
    -c
    8
    -p
    qmc_spindens_dmc_long
    -r
    qmc-ref/qmc_spindens_dmc_long.s001.stat_ref_spindensity.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-rmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_rmc_long.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    20
    -c
    8
    -p
    qmc_spindens_rmc_long
    -r
    qmc-ref/qmc_spindens_rmc_long.s001.stat_ref_spindensity.dat)

  # 1rdm tests

  # no Jastrow tests
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_1rdm_noJ_short.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    1rdm
    -e
    20
    -c
    8
    -p
    qmc_1rdm_noJ_short
    -r
    qmc-ref/qmc_1rdm_noJ_short${OFEXT}.s000.stat_ref_1rdm.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_1rdm_noJ_long.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    1rdm
    -e
    20
    -c
    8
    -p
    qmc_1rdm_noJ_long
    -r
    qmc-ref/qmc_1rdm_noJ_long${OFEXT}.s000.stat_ref_1rdm.dat)

  # two body Jastrow tests
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_1rdm_J2_short.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    1rdm
    -e
    20
    -c
    8
    -p
    qmc_1rdm_J2_short
    -r
    qmc-ref/qmc_1rdm_J2_short${OFEXT}.s000.stat_ref_1rdm.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_1rdm_J2_long.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    1rdm
    -e
    20
    -c
    8
    -p
    qmc_1rdm_J2_long
    -r
    qmc-ref/qmc_1rdm_J2_long${OFEXT}.s000.stat_ref_1rdm.dat)

  # energydensity tests
  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_edens_cell_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    'energydensity,EDcell'
    -e
    20
    -c
    8
    -p
    qmc_edens_cell_short
    -r
    qmc-ref/qmc_edens_cell_short.s000.stat_ref_energydensity.dat)

  simple_run_and_check(
    short-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_edens_vor_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    'energydensity,EDvoronoi'
    -e
    20
    -c
    2
    -p
    qmc_edens_vor_short
    -r
    qmc-ref/qmc_edens_vor_short.s000.stat_ref_energydensity.dat)

  simple_run_and_check(
    short-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_edens_cell_dmc_short${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    'energydensity,EDcell'
    -e
    2
    -c
    8
    -p
    qmc_edens_cell_dmc_short
    -r
    qmc-ref/qmc_edens_cell_dmc_short.s001.stat_ref_energydensity.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_edens_cell_long${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    'energydensity,EDcell'
    -e
    20
    -c
    8
    -p
    qmc_edens_cell_long
    -r
    qmc-ref/qmc_edens_cell_long.s000.stat_ref_energydensity.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_edens_vor_long${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    0
    -q
    'energydensity,EDvoronoi'
    -e
    20
    -c
    2
    -p
    qmc_edens_vor_long
    -r
    qmc-ref/qmc_edens_vor_long.s000.stat_ref_energydensity.dat)

  simple_run_and_check(
    long-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_edens_cell_dmc_long${IFEXT}.in.xml
    ${NMPI}
    ${NOMP}
    check_stats.py
    -s
    1
    -q
    'energydensity,EDcell'
    -e
    2
    -c
    8
    -p
    qmc_edens_cell_dmc_long
    -r
    qmc-ref/qmc_edens_cell_dmc_long.s001.stat_ref_energydensity.dat)

endif()

#
# Deterministic tests
#

# deterministic estimator tests
if(add_estimator_tests AND NOT QMC_MIXED_PRECISION)
  # density tests
  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_dens_short.in.xml
    1
    1
    check_stats.py
    -s
    0
    -q
    density
    -e
    0
    -c
    8
    -p
    det_qmc_dens_short
    -r
    qmc-ref/det_qmc_dens_short.s000.stat_ref_density_1_1.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_dens_short_mw.in.xml
    4
    4
    check_stats.py
    -s
    0
    -q
    density
    -e
    0
    -c
    8
    -p
    det_qmc_dens_short_mw
    -r
    qmc-ref/det_qmc_dens_short_mw.s000.stat_ref_density_4_4.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_dens_short_vmc_dmc.in.xml
    1
    1
    check_stats.py
    -s
    1
    -q
    density
    -e
    1
    -c
    8
    -p
    det_qmc_dens_short_vmc_dmc
    -r
    qmc-ref/det_qmc_dens_short_vmc_dmc.s001.stat_ref_density_1_1${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_dens_short_vmc_dmc.in.xml
    4
    4
    check_stats.py
    -s
    1
    -q
    density
    -e
    1
    -c
    8
    -p
    det_qmc_dens_short_vmc_dmc
    -r
    qmc-ref/det_qmc_dens_short_vmc_dmc.s001.stat_ref_density_4_4${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-rmc-estimator-density
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_dens_short_vmc_rmc.in.xml
    1
    1
    check_stats.py
    -s
    1
    -q
    density
    -e
    1
    -c
    8
    -p
    det_qmc_dens_short_vmc_rmc
    -r
    qmc-ref/det_qmc_dens_short_vmc_rmc.s001.stat_ref_density_1_1.dat)

  # spindensity tests
  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_spindens_short.in.xml
    1
    1
    check_stats.py
    -s
    0
    -q
    spindensity
    -e
    0
    -c
    8
    -p
    det_qmc_spindens_short
    -r
    qmc-ref/det_qmc_spindens_short.s000.stat_ref_spindensity_1_1.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_spindens_short_mw.in.xml
    4
    4
    check_stats.py
    -s
    0
    -q
    spindensity
    -e
    0
    -c
    8
    -p
    det_qmc_spindens_short_mw
    -r
    qmc-ref/det_qmc_spindens_short_mw.s000.stat_ref_spindensity_4_4.dat)

  # SIMPLE_RUN_AND_CHECK(
  #   deterministic-diamondC_1x1x1_pp-vmcbatch-estimator-spindensity
  #   "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  #   det_qmc_spindens_short_vmcbatched.in.xml
  #   1 1
  #   check_stats.py -s 0 -q spindensity -e 0 -c 8 -p det_qmc_spindens_short -r qmc-ref/det_qmc_vmcbatch_spindens_short.s000.stat_ref_spindensity_1_1.dat
  #   )

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_spindens_short_vmc_dmc.in.xml
    1
    1
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    1
    -c
    8
    -p
    det_qmc_spindens_short_vmc_dmc
    -r
    qmc-ref/det_qmc_spindens_short_vmc_dmc.s001.stat_ref_spindensity_1_1${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_spindens_short_vmc_dmc.in.xml
    4
    4
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    1
    -c
    8
    -p
    det_qmc_spindens_short_vmc_dmc
    -r
    qmc-ref/det_qmc_spindens_short_vmc_dmc.s001.stat_ref_spindensity_4_4${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-rmc-estimator-spindensity
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_spindens_short_vmc_rmc.in.xml
    1
    1
    check_stats.py
    -s
    1
    -q
    spindensity
    -e
    1
    -c
    8
    -p
    det_qmc_spindens_short_vmc_rmc
    -r
    qmc-ref/det_qmc_spindens_short_vmc_rmc.s001.stat_ref_spindensity_1_1.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_1rdm_noJ_short.in.xml
    1
    1
    check_stats.py
    -s
    0
    -q
    1rdm
    -e
    0
    -c
    8
    -p
    det_qmc_1rdm_noJ_short
    -r
    qmc-ref/det_qmc_1rdm_noJ_short.s000.stat_ref_1rdm${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_1rdm_J2_short.in.xml
    1
    1
    check_stats.py
    -s
    0
    -q
    1rdm
    -e
    0
    -c
    8
    -p
    det_qmc_1rdm_J2_short
    -r
    qmc-ref/det_qmc_1rdm_J2_short.s000.stat_ref_1rdm${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_1rdm_short_vmc_dmc.in.xml
    1
    1
    check_stats.py
    -s
    1
    -q
    1rdm
    -e
    0
    -c
    8
    -p
    det_qmc_1rdm_short_vmc_dmc
    -r
    qmc-ref/det_qmc_1rdm_short_vmc_dmc.s001.stat_ref_1rdm_1_1${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-1rdm
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_1rdm_short_vmc_dmc.in.xml
    4
    4
    check_stats.py
    -s
    1
    -q
    1rdm
    -e
    0
    -c
    8
    -p
    det_qmc_1rdm_short_vmc_dmc
    -r
    qmc-ref/det_qmc_1rdm_short_vmc_dmc.s001.stat_ref_1rdm_4_4${OFEXT}.dat)

  # momemtum distribution test
  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-estimator-momentum
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_momentum_short.in.xml
    1
    1
    check_stats.py
    -s
    0
    -q
    momentum
    -e
    0
    -c
    5
    -p
    det_qmc_momentum_short
    -r
    qmc-ref/det_qmc_momentum_short.s000.stat_ref_momentum_1_1.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-momentum
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_momentum_short_vmc_dmc.in.xml
    1
    1
    check_stats.py
    -s
    1
    -q
    momentum
    -e
    0
    -c
    5
    -p
    det_qmc_momentum_short_vmc_dmc
    -r
    qmc-ref/det_qmc_momentum_short_vmc_dmc.s001.stat_ref_momentum_1_1${OFEXT}.dat)

  simple_run_and_check(
    deterministic-diamondC_1x1x1_pp-dmc-estimator-momentum
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_momentum_short_vmc_dmc.in.xml
    4
    4
    check_stats.py
    -s
    1
    -q
    momentum
    -e
    0
    -c
    5
    -p
    det_qmc_momentum_short_vmc_dmc
    -r
    qmc-ref/det_qmc_momentum_short_vmc_dmc.s001.stat_ref_momentum_4_4${OFEXT}.dat)

endif()

# deterministic scalar tests

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.65021596 0.00001582")
  list(APPEND DET_DIAMOND_SCALARS "kinetic" "16.03245992 0.00001650")
  list(APPEND DET_DIAMOND_SCALARS "potential" "-26.68267474 0.00000819")
  list(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.64876197 0.00000211")
  list(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "localecp" "-13.45055537 0.00000725")
  list(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "1.19231153 0.00000439")
  list(APPEND DET_DIAMOND_SCALARS "mpc" "-1.45236931 0.00000196")
  list(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_SCALARS "flux" "-14.25663444 0.00005886")
else()
  list(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.65021864 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "kinetic" "16.03245759 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "potential" "-26.68267623 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.64876057 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "localecp" "-13.45055710 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "1.19230887 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "mpc" "-1.45236706 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_SCALARS "flux" "-14.25663001 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short
  det_qmc_short.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_SCALARS # VMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "totenergy" "-10.14834087 0.00001992")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "kinetic" "14.67050590 0.00001980")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "potential" "-24.81885848 0.00003")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "eeenergy" "-3.54361899 0.00001190")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "ionion" "-12.77567050 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "localecp" "-8.35089725 0.00005")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "nonlocalecp" "-0.14867641 0.00002")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "mpc" "-3.24606471 0.00001149")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "flux" "-3.81746292 0.0001")
else()
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "totenergy" "-10.14833448 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "kinetic" "14.67051169 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "potential" "-24.81884682 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "eeenergy" "-3.54361900 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "ionion" "-12.77566767 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "localecp" "-8.35088811 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "nonlocalecp" "-0.14867158 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "mpc" "-3.24606659 0.000001")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_EXCITED_SCALARS "flux" "-3.81745231 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-vmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short_excited
  det_qmc_short_excited.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_EXCITED_SCALARS # VMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.65736954 0.00022771")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "12.60307264 0.00022263")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-23.26040613 0.00001348")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-3.48893504 0.00000280")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-8.52555859 0.00002025")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "1.52975113 0.00000811")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-3.34182284 0.00000307")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "13.85794896 0.00099470")
else()
  list(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.65723898 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "12.60316496 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-23.26040394 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-3.48893469 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-8.52555807 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "1.52975625 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-3.34182165 0.000001")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "13.85780746 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-vmc_sdj-meshf
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short-meshf
  det_qmc_short-meshf.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_MESHF_SCALARS # VMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "totenergy" "-10.43066035 0.00001226")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "kinetic" "11.79557753 0.00001421")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "potential" "-22.22624477 0.00000791")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "eeenergy" "-1.30357070 0.00000629")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "localecp" "-9.17379702 0.00001481")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "nonlocalecp" "1.02680376 0.00000270")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "mpc" "-1.00378374 0.00000754")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "flux" "-10.16222861 0.00006170")
else()
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "totenergy" "-10.43065963 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "kinetic" "11.79557948 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "potential" "-22.22623911 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "eeenergy" "-1.30357077 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "localecp" "-9.17380267 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "nonlocalecp" "1.02680176 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "mpc" "-1.00378421 0.000001")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "samples" "9 0.0")
  list(APPEND DET_DIAMOND_KSPACE_SCALARS "flux" "-10.16222139 0.000001")
endif()

if(NOT ENABLE_OFFLOAD)
  qmc_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc_sdj_kspace
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_short_kspace
    det_qmc_short_kspace.in.xml
    1
    1
    TRUE
    0
    DET_DIAMOND_KSPACE_SCALARS # VMC
  )
endif()

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.52331645 0.00001129")
  list(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "4.49057680 0.00001257")
  list(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-15.01388990 0.00000253")
  list(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.83205956 0.00000157")
  list(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.00159178 0.00000147")
  list(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "0.93131260 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.40050428 0.00000225")
  list(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
else()
  list(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.52330287 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "4.49059156 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-15.01389443 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.83205797 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.00159889 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "0.93131392 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.40050042 0.000001")
  list(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short_vmc_dmc
  det_qmc_short_vmc_dmc.in.xml
  1
  1
  TRUE
  1
  DET_DIAMOND_DMC_SCALARS # DMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-9.83668425 0.00001562")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "kinetic" "9.63546916 0.00001143")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "potential" "-19.47214772 0.00000463")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-2.70238228  0.00000125")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "ionion" "-12.77567050 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "localecp" "-4.36957276 0.0000047")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "0.37547763 0.000001 ")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "mpc" "-2.42116567 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "samples" "9 0.0")
else()
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-9.83667502 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "kinetic" "9.63547298 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "potential" "-19.47214800 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-2.70238234 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "ionion" "-12.77566767 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "localecp" "-4.36957578 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "0.37547753 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "mpc" "-2.421165687 0.000001")
  list(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "samples" "9 0.0")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-dmc_sdj_excited
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short_vmc_dmc_excited
  det_qmc_short_vmc_dmc_excited.in.xml
  1
  1
  TRUE
  1
  DET_DIAMOND_DMC_EXCITED_SCALARS # DMC
)

if(QMC_MIXED_PRECISION)
  # DMC cycle 1
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "totenergy" "-9.53121536 0.0000139")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "kinetic" "9.00930034 0.00001647")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "potential" "-18.54051114 0.00000446")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "eeenergy" "-2.11983595 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "localecp" "-4.36659965 0.00000448")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "nonlocalecp" "0.72159739 0.00000178")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "mpc" "-1.93675985 0.00000122")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "samples" "9 0.0")
  # DMC cycle 2
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "totenergy" "-10.56881667 0.000027")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "kinetic" "14.41362921 0.00003329")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "potential" "-24.98245366 0.00003051")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "eeenergy" "-2.89879032 0.00000133")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "localecp" "-15.41387689 0.00001471")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "nonlocalecp" "6.10588419 0.00001448")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "mpc" "-2.61905734 0.00000179")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "samples" "9 0.0")
  # DMC cycle 3
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "totenergy" "-8.93222965 0.00001861")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "kinetic" "19.70827934 0.00003064")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "potential" "-28.64051062 0.00001851")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "eeenergy" "-1.15369571 0.00000776")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "localecp" "-16.18244873 0.00001536")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "nonlocalecp" "1.47130834 0.00000689")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "mpc" "-0.80671517 0.00000857")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "samples" "9 0.0")
else()
  # DMC cycle 1
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "totenergy" "-9.53119806 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "kinetic" "9.00930862 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "potential" "-18.54050668 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "eeenergy" "-2.11983454 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "localecp" "-4.36660098 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "nonlocalecp" "0.72159631 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "mpc" "-1.93675646 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "samples" "9 0.0")
  # DMC cycle 2
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "totenergy" "-10.56879581 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "kinetic" "14.41364745 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "potential" "-24.98244325 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "eeenergy" "-2.89879008 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "localecp" "-15.41389327 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "nonlocalecp" "6.10590756 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "mpc" "-2.61905530 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "samples" "9 0.0")
  # DMC cycle 3
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "totenergy" "-8.93222899 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "kinetic" "19.70828644 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "potential" "-28.64051544 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "eeenergy" "-1.15369726 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "localecp" "-16.18246211 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "nonlocalecp" "1.47131139 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "mpc" "-0.80671598 0.000001")
  list(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "samples" "9 0.0")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-multi_dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short_vmc_multi_dmc
  det_qmc_short_vmc_multi_dmc.in.xml
  1
  1
  TRUE
  1
  DET_DIAMOND_DMC_MULTI1_SCALARS # DMC cycle 1
  2
  DET_DIAMOND_DMC_MULTI2_SCALARS # DMC cycle 2
  3
  DET_DIAMOND_DMC_MULTI3_SCALARS # DMC cycle 3
)

if(QMC_MIXED_PRECISION)
  # VMC cycle 1
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "totenergy" "-10.51887495 0.00000925")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "kinetic" "12.21076917 0.00000858")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "potential" "-22.72964386 0.00000257")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "eeenergy" "-1.39227245 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "localecp" "-8.27168446 0.0000035")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "nonlocalecp" "-0.29001602 0.00000126")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "mpc" "-0.95019080 0.00000112")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # VMC cycle 2
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "totenergy" "-10.57006353 0.00001796")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "kinetic" "13.81347725 0.0000143")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "potential" "-24.38354103 0.00000366")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "eeenergy" "-2.08653505 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "ionion" "-12.77567000  0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "localecp" "-10.59945535 0.00000897")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "nonlocalecp" "1.07811994 0.00000743")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "mpc" "-2.00098290 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "samples" "18 0.0")
  # DMC cycle 1
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "totenergy" "-10.53543619 0.00000979")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "kinetic" "7.04049333 0.00000959")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "potential" "-17.57593026 0.00000443")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "eeenergy" "-1.15002939 0.00000344")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "localecp" "-3.78048184 0.00000184")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "nonlocalecp" "0.13025148 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "mpc" "-1.01499425 0.00000365")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # DMC cycle2
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "totenergy" "-9.72346878 0.0000083")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "kinetic" "11.45079701 0.00000847")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "potential" "-21.17426291 0.00000232")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "eeenergy" "-1.73937136 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "localecp" "-8.33546910 0.00000206")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "nonlocalecp" "1.67624709 0.00000168")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "mpc" "-1.50889162 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "samples" "18 0.0")
else()
  # VMC cycle 1
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "totenergy" "-10.51886527 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "kinetic" "12.21077547 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "potential" "-22.72964075 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "eeenergy" "-1.39227187 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "localecp" "-8.27169251 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "nonlocalecp" "-0.29000890 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "mpc" "-0.95018933 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # VMC cycle 2
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "totenergy" "-10.57006352 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "kinetic" "13.81348242 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "potential" "-24.38354594 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "eeenergy" "-2.08653334 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "localecp" "-10.59945467 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "nonlocalecp" "1.07810954 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "mpc" "-2.00098039 0.000001")
  list(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "samples" "18 0.0")
  # DMC cycle 1
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "totenergy" "-10.53543619 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "kinetic" "7.04048475 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "potential" "-17.57592094 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "eeenergy" "-1.15002397 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "localecp" "-3.78048098 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "nonlocalecp" "0.13025148 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "mpc" "-1.01498843 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # DMC cycle 2
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "totenergy" "-9.72346093 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "kinetic" "11.45079979 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "potential" "-21.17426072 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "eeenergy" "-1.73936902 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "localecp" "-8.33547227 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "nonlocalecp" "1.67624804 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "mpc" "-1.50888780 0.000001")
  list(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "samples" "18 0.0")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-vmc-dmc-mwalkers_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short_vmc_dmc_mwalkers
  det_qmc_short_vmc_dmc_mwalkers.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_VMC_MWALKERS1_SCALARS # VMC cycle 1
  1
  DET_DIAMOND_VMC_MWALKERS2_SCALARS # VMC cycle 2
  2
  DET_DIAMOND_DMC_MWALKERS1_SCALARS # DMC cycle 1
  3
  DET_DIAMOND_DMC_MWALKERS2_SCALARS # DMC cycle 2
)

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-vmc-dmc-mwalkers-checks_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_short_vmc_dmc_mwalkers
  det_qmc_short_vmc_dmc_mwalkers_with_checks.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_VMC_MWALKERS1_SCALARS # VMC cycle 1
  1
  DET_DIAMOND_VMC_MWALKERS2_SCALARS # VMC cycle 2
  2
  DET_DIAMOND_DMC_MWALKERS1_SCALARS # DMC cycle 1
  3
  DET_DIAMOND_DMC_MWALKERS2_SCALARS # DMC cycle 2
)

if((NOT QMC_MIXED_PRECISION) AND (NOT QMC_COMPLEX))
  # DMC cycle 1
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "totenergy" "-10.77242244 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "kinetic" "11.31489895 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "potential" "-22.08732139 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "eeenergy" "-3.20777721 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "localecp" "-7.31306387 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "nonlocalecp" "1.20918714 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "mpc" "-3.10215768 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "samples" "15 0.0")
  # DMC cycle 2
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "totenergy" "-10.18494655 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "kinetic" "13.02680354 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "potential" "-23.21175009 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "eeenergy" "-2.77542078 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "localecp" "-8.59996640 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "nonlocalecp" "0.93930455 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "mpc" "-2.48119177 0.000001")
  list(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "samples" "18 0.0")

  qmc_run_and_check(
    deterministic-diamondC_1x1x1_pp-vmc-dmc-multi_ts_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    det_qmc_short_vmc_dmc_lg_ts
    det_qmc_short_vmc_dmc_lg_ts.in.xml
    1
    1
    TRUE
    1
    DET_DIAMOND_DMC_LGTS1_SCALARS # DMC cycle 1
    2
    DET_DIAMOND_DMC_LGTS2_SCALARS # DMC cycle 2
  )
endif()

if(QMC_MIXED_PRECISION)
  #v0
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "totenergy" "-10.44253935 0.00001861")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "kinetic" "11.35705140 0.00003064")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "potential" "-21.79959076 0.00000232")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "eeenergy" "-0.69361370 0.00000776")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "localecp" "-9.65190864 0.00001536")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "nonlocalecp" "1.32159905 0.00000689")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "mpc" "-0.23519803 0.00000857")
  #v1
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "totenergy" "-10.49290122 0.00001861")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "kinetic" "12.78036761 0.00003064")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "potential" "-23.27326883 0.00000232")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "eeenergy" "-1.36300017 0.00000776")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "localecp" "-10.21261249 0.00001536")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "nonlocalecp" "1.07801130 0.00000689")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "mpc" "-0.98142311 0.00000857")
  #v3
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "totenergy" "-10.51064257 0.00001861")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "kinetic" "14.08671956 0.00003064")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "potential" "-24.59736212 0.00000232")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "eeenergy" "-1.90175515 0.00000776")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "localecp" "-11.18722424 0.00001536")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "nonlocalecp" "1.26728473 0.00000689")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "mpc" "-1.59078038 0.00000857")
else()
  #v0
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "totenergy" "-10.44253935 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "kinetic" "11.35705140 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "potential" "-21.79959076 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "eeenergy" "-0.69361370 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "localecp" "-9.65190864 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "nonlocalecp" "1.32159905 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_SCALARS "mpc" "-0.23519803 0.000001")
  #v1
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "totenergy" "-10.49290122 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "kinetic" "12.78036761 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "potential" "-23.27326883 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "eeenergy" "-1.36300017 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "localecp" "-10.21261249 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "nonlocalecp" "1.07801130 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_SCALARS "mpc" "-0.98142311 0.000001")
  #v3
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "totenergy" "-10.51064257 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "kinetic" "14.08671956 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "potential" "-24.59736212 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "eeenergy" "-1.90175515 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "localecp" "-11.18722424 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "nonlocalecp" "1.26728473 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_SCALARS "mpc" "-1.59078038 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-dmc-tmoves_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_vmc_dmc_tm
  det_qmc_vmc_dmc_tm.in.xml
  1
  1
  TRUE
  1
  DET_DIAMOND_DMC_TM1_SCALARS # DMC v0
  2
  DET_DIAMOND_DMC_TM2_SCALARS # DMC v1
  3
  DET_DIAMOND_DMC_TM3_SCALARS # DMC v3
)

if(QMC_MIXED_PRECISION)
  #v0
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "totenergy" "-10.34857080 0.00001861")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "kinetic" "9.94328567 0.00003064")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "potential" "-20.29186355 0.00000232")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "eeenergy" "-2.70623197 0.00000776")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "localecp" "-5.26283000 0.00001536")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "nonlocalecp" "0.45287649 0.00000689")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "mpc" "-2.43479306 0.00000857")
  #v1
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "totenergy" "-10.47605841 0.00001861")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "kinetic" "10.81438943 0.00003064")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "potential" "-21.29044786 0.00000232")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "eeenergy" "-2.90631955 0.00000776")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "localecp" "-5.56045812 0.00001536")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "nonlocalecp" "-0.04800030 0.00000689")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "mpc" "-2.64915619 0.00000857")
  #v3
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "totenergy" "-10.29166891 0.00001861")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "kinetic" "12.08875246 0.00003064")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "potential" "-22.38042133 0.00000232")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "eeenergy" "-2.82248174 0.00000776")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "localecp" "-6.53921012 0.00001536")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "nonlocalecp" "-0.24305956 0.00000689")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "mpc" "-2.57378019 0.00000857")
else()
  #v0
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "totenergy" "-10.74426355 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "kinetic" "9.36905833 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "potential" "-20.11332161 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "eeenergy" "-2.55746488 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "localecp" "-5.08312451 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "nonlocalecp" "0.30293477 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM1_BATCH_SCALARS "mpc" "-2.31725901 0.000001")
  #v1
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "totenergy" "-10.78859676 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "kinetic" "11.17581012 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "potential" "-21.96440688 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "eeenergy" "-2.55624212 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "localecp" "-6.99272066 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "nonlocalecp" "0.36022260 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM2_BATCH_SCALARS "mpc" "-2.31923312 0.000001")
  #v3
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "totenergy" "-10.77390772 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "kinetic" "13.76835233 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "potential" "-24.54226005 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "eeenergy" "-2.69128089 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "ionion" "-12.77567000 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "localecp" "-9.09412623 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "nonlocalecp" "0.01881453 0.000001")
  list(APPEND DET_DIAMOND_DMC_TM3_BATCH_SCALARS "mpc" "-2.46961810 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-dmcbatch-tmoves_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_vmcbatch_dmcbatch_tm
  det_qmc_vmcbatch_dmcbatch_tm.in.xml
  1
  1
  TRUE
  1
  DET_DIAMOND_DMC_TM1_BATCH_SCALARS # DMC v0
  2
  DET_DIAMOND_DMC_TM2_BATCH_SCALARS # DMC v1
  3
  DET_DIAMOND_DMC_TM3_BATCH_SCALARS # DMC v3
)

# traces tests
check_python_reqs("numpy;h5py" diamond-traces add_traces_tests)
if(add_traces_tests)
  if(QMC_MIXED_PRECISION)
    # traces with only scalars on
    simple_run_and_check(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars
      "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars.in.xml
      1
      1
      check_traces.py
      -p
      qmc_trace_scalars
      -s
      '0,1'
      -m
      'vmc,dmc'
      --dmc_steps_exclude=1
      --tol=2e-6)
  else()
    simple_run_and_check(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars
      "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars.in.xml
      1
      1
      check_traces.py
      -p
      qmc_trace_scalars
      -s
      '0,1'
      -m
      'vmc,dmc'
      --dmc_steps_exclude=1)
  endif()

  if(QMC_MIXED_PRECISION)
    # traces with only selective scalars
    simple_run_and_check(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars_sel
      "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars_selective.in.xml
      1
      1
      check_traces.py
      -p
      qmc_trace_scalars_selective
      -s
      '0,1'
      -m
      'vmc,dmc'
      -q
      'Kinetic,ElecElec'
      --dmc_steps_exclude=1
      --tol=1e-6)
  else()
    simple_run_and_check(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars_sel
      "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars_selective.in.xml
      1
      1
      check_traces.py
      -p
      qmc_trace_scalars_selective
      -s
      '0,1'
      -m
      'vmc,dmc'
      -q
      'Kinetic,ElecElec'
      --dmc_steps_exclude=1)
  endif()
endif()

if(NOT QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_PARAM uu_0 -1.010849732e-05 0.00001)
  list(APPEND DET_DIAMOND_PARAM uu_1 1.221957505 0.00010)
  list(APPEND DET_DIAMOND_PARAM ud_0 1.094114467e-07 0.00001)
  list(APPEND DET_DIAMOND_PARAM ud_1 0.0001130922027 0.00001)
  list(APPEND DET_DIAMOND_PARAM eC_2 0.0 0.00001)
  list(APPEND DET_DIAMOND_PARAM eC_3 8.196343356 0.00001)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    deterministic-diamondC_1x1x1_pp-param-grad
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    PREFIX
    det_qmc_short_param_grad.param
    INPUT_FILE
    det_qmc_param_grad.xml
    PROCS
    1
    THREADS
    1
    SERIES
    0
    EQUILIBRATION
    0
    SCALAR_VALUES
    DET_DIAMOND_PARAM)

  list(APPEND DET_DIAMOND_PARAM_LEGACY uu_0 1.1546e-08 0.00001)
  list(APPEND DET_DIAMOND_PARAM_LEGACY uu_1 5.9577e-06 0.00010)
  list(APPEND DET_DIAMOND_PARAM_LEGACY ud_0 0.00000051 0.00001)
  list(APPEND DET_DIAMOND_PARAM_LEGACY ud_1 -0.00568847 0.00001)
  list(APPEND DET_DIAMOND_PARAM_LEGACY eC_2 1.8608 0.00002)
  list(APPEND DET_DIAMOND_PARAM_LEGACY eC_3 4.7107 0.00003)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    deterministic-diamondC_1x1x1_pp-param-grad-legacy-driver
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    PREFIX
    det_qmc_short_param_grad.param
    INPUT_FILE
    det_qmc_param_grad_legacy_driver.xml
    PROCS
    1
    THREADS
    1
    SERIES
    0
    EQUILIBRATION
    0
    SCALAR_VALUES
    DET_DIAMOND_PARAM_LEGACY)
endif()

#
# Basic input file checking tests: malformed, missing input, missing files
#
qmc_run_and_check(
  deterministic-input_check_broken_xml
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_broken_xml
  det_qmc_input_check_broken_xml.in.xml
  1
  1
  FALSE)

qmc_run_and_check(
  deterministic-input_check_no_ham
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_no_ham
  det_qmc_input_check_no_ham.in.xml
  1
  1
  FALSE)

qmc_run_and_check(
  deterministic-input_check_no_pp
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_no_pp
  det_qmc_input_check_no_pp.in.xml
  1
  1
  FALSE)

qmc_run_and_check(
  deterministic-input_check_no_pset
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_no_pset
  det_qmc_input_check_no_pset.in.xml
  1
  1
  FALSE)

qmc_run_and_check(
  deterministic-input_check_no_qmc
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_no_qmc
  det_qmc_input_check_no_qmc.in.xml
  1
  1
  TRUE)

qmc_run_and_check(
  deterministic-input_check_no_splineh5
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_no_splineh5
  det_qmc_input_check_no_splineh5.in.xml
  1
  1
  FALSE)

qmc_run_and_check(
  deterministic-input_check_no_wf
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_no_wf
  det_qmc_input_check_no_wf.in.xml
  1
  1
  FALSE)

qmc_run_and_check(
  deterministic-input_check_not_xml
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
  det_qmc_input_check_not_xml
  det_qmc_input_check_not_xml.in.xml
  1
  1
  FALSE)
