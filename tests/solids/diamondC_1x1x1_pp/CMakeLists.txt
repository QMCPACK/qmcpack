
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.491445 0.002")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.434392 0.018")
  LIST(APPEND DIAMOND_SCALARS "potential" "-21.925822 0.018")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.700534 0.0043")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566752 0.001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.046740 0.020")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.597119 0.0056")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
  LIST(APPEND DIAMOND_SCALARS "mpc" "-2.453044 0.004431")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_hybridrep_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_hybridrep.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short-meshf
                    qmc_short-meshf.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

IF(NOT QMC_CUDA AND NOT ENABLE_OFFLOAD)
  LIST(APPEND DIAMOND_KSPACE_SCALARS "totenergy" "-10.500719 0.001769")
  LIST(APPEND DIAMOND_KSPACE_SCALARS "variance"  "0.312264 0.028662")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj_kspace
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_kspace
                    qmc_short_kspace.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_KSPACE_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj_kspace
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_kspace_4_4
                    qmc_short_kspace_4_4.in.xml
                    4 4
                    TRUE
                    0 DIAMOND_KSPACE_SCALARS # VMC
                    )
ELSE()
  MESSAGE_VERBOSE("Skipping k-space Jastrow tests because they are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()


# Reference OPT run in qmc-ref
  LIST(APPEND DIAMOND_OPT_SCALARS "totenergy"   "-10.49370 0.0024")
  LIST(APPEND DIAMOND_OPT_SCALARS "kinetic"     "11.5504 0.021")
  LIST(APPEND DIAMOND_OPT_SCALARS "localecp"    "-7.1979 0.023")
  LIST(APPEND DIAMOND_OPT_SCALARS "nonlocalecp" "0.61540 0.0068")
  LIST(APPEND DIAMOND_OPT_SCALARS "samples"     "128000 0.0")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-opt_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_opt.in.xml
                    1 16
                    ${MP_SUCCESS}
                    3 DIAMOND_OPT_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND DIAMOND_DMC_SCALARS "totenergy"   "-10.5316 0.0024")
  LIST(APPEND DIAMOND_DMC_SCALARS "kinetic"      "11.4857 0.0231")
  LIST(APPEND DIAMOND_DMC_SCALARS "potential"   "-22.0170 0.0239")
  LIST(APPEND DIAMOND_DMC_SCALARS "localecp"     "-7.1518 0.0323")
  LIST(APPEND DIAMOND_DMC_SCALARS "nonlocalecp"  "0.62688 0.0080")
  LIST(APPEND DIAMOND_DMC_SCALARS "eeenergy"    "-2.71641 0.0073")
  LIST(APPEND DIAMOND_DMC_SCALARS "mpc"    "-2.47994 0.00779")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND_DMC_SCALARS # DMC
                    )

# Excited state
  LIST(APPEND DIAMOND_EXCITED_SCALARS "totenergy"  "-9.997214 0.002156")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "kinetic"    "11.454674 0.017147")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "potential" "-21.451888 0.017343")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "eeenergy"   "-2.790777 0.004311")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "ionion"    "-12.775668 0.00001")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "localecp"   "-6.493792 0.019499")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "nonlocalecp" "0.608349 0.005683")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "samples" "128000 0.0")
  LIST(APPEND DIAMOND_EXCITED_SCALARS "mpc" "-2.504948 0.004507")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_excited
                    qmc_short_excited.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_EXCITED_SCALARS # VMC
                    )

  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "totenergy"   "-10.057283 0.003216")
  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "kinetic"      "11.598366 0.027848")
  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "potential"   "-21.655650 0.029074")
  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "localecp"     "-6.711209 0.035888")
  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp"   "0.647451 0.008633")
  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "eeenergy"     "-2.816225 0.007467")
  LIST(APPEND DIAMOND_DMC_EXCITED_SCALARS "mpc"    "-2.548944 0.00801")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-dmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_vmc_dmc_excited
                    qmc_short_vmc_dmc_excited.in.xml
                    1 16
                    TRUE
                    1 DIAMOND_DMC_EXCITED_SCALARS # DMC
                    )

IF(NOT QMC_CUDA)
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc-dmc-allp_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_vmc_dmc_allp
                    qmc_short_vmc_dmc_allp.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC without drift
                    1 DIAMOND_SCALARS # VMC with drift
                    2 DIAMOND_DMC_SCALARS # DMC
                    )
ELSE()
  MESSAGE("Skipping diamondC_1x1x1_pp all particle move VMC/DMC tests because they are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

# Coverage test - shorter version of dmc test
  COVERAGE_RUN("diamondC_1x1x1_pp-vmc_sdj-1-1"
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              qmc_short_vmc_dmc_coverage.in.xml)

# Test time limit
  CPU_LIMIT_RUN("diamondC_1x1x1_pp-vmc-1-1"
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              180
              qmc_cpu_limit_vmc.xml)

  CPU_LIMIT_RUN("diamondC_1x1x1_pp-vmc-4-4"
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              4 4
              120
              qmc_cpu_limit_vmc.xml)

  CPU_LIMIT_RUN("diamondC_1x1x1_pp-dmc-1-1"
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              240
              qmc_cpu_limit_dmc.xml)

  CPU_LIMIT_RUN("diamondC_1x1x1_pp-dmc-4-4"
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              4 4
              120
              qmc_cpu_limit_dmc.xml)


#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.491445 0.000065")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "122880000 0.0")
  LIST(APPEND LONG_DIAMOND_SCALARS "flux" "0.0 0.03")
  LIST(APPEND LONG_DIAMOND_SCALARS "mpc" "-2.453044 0.000143")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long-meshf
                    qmc_long-meshf.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

IF(NOT QMC_CUDA)
  # VMC ref
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "totenergy" "-10.491445 0.0065")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "kinetic" "11.434392 0.006")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "potential" "-21.925822 0.006")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "eeenergy" "-2.700534 0.0015")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "ionion" "-12.77566752 0.0004")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "localecp" "-7.046740 0.0065")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "nonlocalecp" "0.597119 0.0019")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "mpc" "-2.453044 0.0014")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "samples" "1280000 0.0")
  #DMC ref
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "totenergy"   "-10.5316 0.0008")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "kinetic"      "11.4857 0.0077")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "potential"   "-22.0170 0.0080")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "localecp"     "-7.1518 0.0108")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "nonlocalecp"  "0.62688 0.0027")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "eeenergy"    "-2.71641 0.0028")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "mpc" "-2.47943 0.00268")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc-dmc-allp_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_vmc_dmc_allp
                    qmc_long_vmc_dmc_allp.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_ALLP_SCALARS # VMC without drift
                    1 LONG_DIAMOND_ALLP_SCALARS # VMC with drift
                    2 LONG_DIAMOND_DMC_ALLP_SCALARS # DMC
                    )
ENDIF()

IF(NOT QMC_CUDA AND NOT ENABLE_OFFLOAD)
  LIST(APPEND LONG_DIAMOND_KSPACE_SCALARS "totenergy" "-10.500719 0.001769")
  LIST(APPEND LONG_DIAMOND_KSPACE_SCALARS "variance"  "0.312264 0.028662")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj_kspace
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_kspace
                    qmc_long_kspace.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_KSPACE_SCALARS # VMC
                    )
ELSE()
  MESSAGE_VERBOSE("Skipping k-space Jastrow tests because they are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.531583 0.000815")
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "mpc" "-2.47994 0.00246")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

# Excited state
  LIST(APPEND LONG_DIAMOND_EXCITED_SCALARS "totenergy" "-9.997214 0.000073")
  LIST(APPEND LONG_DIAMOND_EXCITED_SCALARS "samples" "122880000 0.0")
  LIST(APPEND LONG_DIAMOND_EXCITED_SCALARS "flux" "-0.006826 0.040699")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_excited
                    qmc_long_excited.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_EXCITED_SCALARS # VMC
                    )

  LIST(APPEND LONG_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-10.057283 0.001061")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_vmc_dmc_excited
                    qmc_long_vmc_dmc_excited.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_EXCITED_SCALARS # DMC
                    )

# estimator tests
INCLUDE("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
set(stat_python_reqs numpy;h5py)
CHECK_PYTHON_REQS(stat_python_reqs diamond-estimator add_estimator_tests)
if (NOT QMC_CUDA)
  set(NMPI 4)
  set(NOMP 4)
  set(IFEXT "")
else()
  set(NMPI 1)
  set(NOMP 16)
  set(IFEXT "_gpu")
endif()
if (NOT QMC_COMPLEX)
  set(OFEXT "_real")
else()
  set(OFEXT "_comp")
endif()
if (add_estimator_tests)
  # density tests
  SIMPLE_RUN_AND_CHECK(
    short-diamondC_1x1x1_pp-vmc-estimator-density
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_short${IFEXT}.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 0 -q density -e 20 -c 8 -p qmc_dens_short -r qmc-ref/qmc_dens_short.s000.stat_ref_density.dat
    )

  SIMPLE_RUN_AND_CHECK(
    short-diamondC_1x1x1_pp-dmc-estimator-density
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_dmc_short.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 1 -q density -e 20 -c 8 -p qmc_dens_dmc_short -r qmc-ref/qmc_dens_dmc_short.s001.stat_ref_density.dat
    )

  SIMPLE_RUN_AND_CHECK(
    long-diamondC_1x1x1_pp-vmc-estimator-density
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_long${IFEXT}.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 0 -q density -e 20 -c 8 -p qmc_dens_long -r qmc-ref/qmc_dens_long.s000.stat_ref_density.dat
    )

  SIMPLE_RUN_AND_CHECK(
    long-diamondC_1x1x1_pp-dmc-estimator-density
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_dens_dmc_long.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 1 -q density -e 20 -c 8 -p qmc_dens_dmc_long -r qmc-ref/qmc_dens_dmc_long.s001.stat_ref_density.dat
    )


  # spindensity tests
  SIMPLE_RUN_AND_CHECK(
    short-diamondC_1x1x1_pp-vmc-estimator-spindensity
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_short${IFEXT}.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 0 -q spindensity -e 20 -c 8 -p qmc_spindens_short -r qmc-ref/qmc_spindens_short.s000.stat_ref_spindensity.dat
    )

  SIMPLE_RUN_AND_CHECK(
    short-diamondC_1x1x1_pp-dmc-estimator-spindensity
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_dmc_short.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 1 -q spindensity -e 2 -c 8 -p qmc_spindens_dmc_short -r qmc-ref/qmc_spindens_dmc_short.s001.stat_ref_spindensity.dat
    )

  SIMPLE_RUN_AND_CHECK(
    long-diamondC_1x1x1_pp-vmc-estimator-spindensity
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_long${IFEXT}.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 0 -q spindensity -e 20 -c 8 -p qmc_spindens_long -r qmc-ref/qmc_spindens_long.s000.stat_ref_spindensity.dat
    )

  SIMPLE_RUN_AND_CHECK(
    long-diamondC_1x1x1_pp-dmc-estimator-spindensity
    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
    qmc_spindens_dmc_long.in.xml
    ${NMPI} ${NOMP}
    check_stats.py -s 1 -q spindensity -e 20 -c 8 -p qmc_spindens_dmc_long -r qmc-ref/qmc_spindens_dmc_long.s001.stat_ref_spindensity.dat
    )

  # 1rdm tests
  if (NOT QMC_CUDA)
    # no Jastrow tests
    SIMPLE_RUN_AND_CHECK(
      short-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_1rdm_noJ_short.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 1rdm -e 20 -c 8 -p qmc_1rdm_noJ_short -r qmc-ref/qmc_1rdm_noJ_short${OFEXT}.s000.stat_ref_1rdm.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-diamondC_1x1x1_pp-vmc-noJ-estimator-1rdm
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_1rdm_noJ_long.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 1rdm -e 20 -c 8 -p qmc_1rdm_noJ_long -r qmc-ref/qmc_1rdm_noJ_long${OFEXT}.s000.stat_ref_1rdm.dat
      )

    # two body Jastrow tests
    SIMPLE_RUN_AND_CHECK(
      short-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_1rdm_J2_short.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 1rdm -e 20 -c 8 -p qmc_1rdm_J2_short -r qmc-ref/qmc_1rdm_J2_short${OFEXT}.s000.stat_ref_1rdm.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-diamondC_1x1x1_pp-vmc-J2-estimator-1rdm
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_1rdm_J2_long.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 1rdm -e 20 -c 8 -p qmc_1rdm_J2_long -r qmc-ref/qmc_1rdm_J2_long${OFEXT}.s000.stat_ref_1rdm.dat
      )
  else()
    message("Skipping diamond 1rdm test because CUDA is not supported")
  endif()

  if (NOT QMC_CUDA)
    # energydensity tests
    SIMPLE_RUN_AND_CHECK(
      short-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_edens_cell_short${IFEXT}.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 'energydensity,EDcell' -e 20 -c 8 -p qmc_edens_cell_short -r qmc-ref/qmc_edens_cell_short.s000.stat_ref_energydensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      short-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_edens_vor_short${IFEXT}.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 'energydensity,EDvoronoi' -e 20 -c 2 -p qmc_edens_vor_short -r qmc-ref/qmc_edens_vor_short.s000.stat_ref_energydensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      short-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_edens_cell_dmc_short${IFEXT}.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 1 -q 'energydensity,EDcell' -e 2 -c 8 -p qmc_edens_cell_dmc_short -r qmc-ref/qmc_edens_cell_dmc_short.s001.stat_ref_energydensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-diamondC_1x1x1_pp-vmc-estimator-energydensity-cell
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_edens_cell_long${IFEXT}.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 'energydensity,EDcell' -e 20 -c 8 -p qmc_edens_cell_long -r qmc-ref/qmc_edens_cell_long.s000.stat_ref_energydensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-diamondC_1x1x1_pp-vmc-estimator-energydensity-voronoi
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_edens_vor_long${IFEXT}.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 0 -q 'energydensity,EDvoronoi' -e 20 -c 2 -p qmc_edens_vor_long -r qmc-ref/qmc_edens_vor_long.s000.stat_ref_energydensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-diamondC_1x1x1_pp-dmc-estimator-energydensity-cell
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_edens_cell_dmc_long${IFEXT}.in.xml
      ${NMPI} ${NOMP}
      check_stats.py -s 1 -q 'energydensity,EDcell' -e 2 -c 8 -p qmc_edens_cell_dmc_long -r qmc-ref/qmc_edens_cell_dmc_long.s001.stat_ref_energydensity.dat
      )

  else()
    message("Skipping diamond 1rdm test because CUDA is not supported")
  endif()


endif()

#
# Deterministic tests
#
IF(NOT QMC_CUDA)

  # deterministic estimator tests
  if (add_estimator_tests AND NOT QMC_MIXED_PRECISION)
    # density tests
    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-estimator-density
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      det_qmc_dens_short.in.xml
      1 1
      check_stats.py -s 0 -q density -e 0 -c 8 -p det_qmc_dens_short -r qmc-ref/det_qmc_dens_short.s000.stat_ref_density_1_1.dat
      )

    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-estimator-density
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      det_qmc_dens_short_mw.in.xml
      4 4
      check_stats.py -s 0 -q density -e 0 -c 8 -p det_qmc_dens_short_mw -r qmc-ref/det_qmc_dens_short_mw.s000.stat_ref_density_4_4.dat
      )
  
    # spindensity tests
    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-estimator-spindensity
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      det_qmc_spindens_short.in.xml
      1 1
      check_stats.py -s 0 -q spindensity -e 0 -c 8 -p det_qmc_spindens_short -r qmc-ref/det_qmc_spindens_short.s000.stat_ref_spindensity_1_1.dat
      )

    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-estimator-spindensity
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      det_qmc_spindens_short_mw.in.xml
      4 4
      check_stats.py -s 0 -q spindensity -e 0 -c 8 -p det_qmc_spindens_short_mw -r qmc-ref/det_qmc_spindens_short_mw.s000.stat_ref_spindensity_4_4.dat
      )
  endif()


  # deterministic scalar tests

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.65155712 0.00001582")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "16.03245992 0.00001650")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-26.68401704 0.00000819")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.64876197 0.00000211")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77567016 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-13.45055537 0.00000725")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "1.19097047 0.00000439")
   LIST(APPEND DET_DIAMOND_SCALARS "mpc" "-1.45236931 0.00000196")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_SCALARS "flux" "-14.25663444 0.00005886")
 ELSE()
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.65156025 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "16.03245759 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-26.68401784 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.64876057 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566756 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-13.45055710 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "1.19096740 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "mpc" "-1.45236706 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_SCALARS "flux" "-14.25663001 0.000001")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short
                    det_qmc_short.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_SCALARS # VMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "totenergy" "-9.41184528 0.00001992")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "kinetic" "12.45900142 0.00001980")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "potential" "-21.87084670 0.00001931")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "eeenergy" "-2.52408721 0.00001190")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "ionion" "-12.77567016 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "localecp" "-6.28770825 0.00001516")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "nonlocalecp" "-0.28338106 0.00000595")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "mpc" "-2.25993608 0.00001149")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "flux" "-14.97659183 0.00004665")
 ELSE()
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "totenergy" "-9.41183124 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "kinetic" "12.45901380 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "potential" "-21.87084504 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "eeenergy" "-2.52408475 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "ionion" "-12.77566744 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "localecp" "-6.28771209 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "nonlocalecp" "-0.28338076 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "mpc" "-2.25993371 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "flux" "-14.97660158 0.000001")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_excited
                    det_qmc_short_excited.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_EXCITED_SCALARS # VMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.65627080 0.00022771")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "12.60307264 0.00022263")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-23.25934345 0.00001348")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-3.48893504 0.00000280")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77567016 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-8.52555859 0.00002025")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "1.53082035 0.00000811")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-3.34182284 0.00000307")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "13.85794896 0.00099470")
 ELSE()
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.65617656 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "12.60316496 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-23.25934152 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-3.48893469 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77566756 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-8.52555807 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "1.53081881 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-3.34182165 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "13.85780746 0.000001")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short-meshf
                    det_qmc_short-meshf.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_MESHF_SCALARS # VMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "totenergy" "-10.42800920 0.00001226")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "kinetic" "11.79557753 0.00001421")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "potential" "-22.22358673 0.00000791")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "eeenergy" "-1.30357070 0.00000629")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "ionion" "-12.77567016 0.000001")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "localecp" "-9.17379702 0.00001481")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "nonlocalecp" "1.02945115 0.00000270")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "mpc" "-1.00378374 0.00000754")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "flux" "-10.16222861 0.00006170")
 ELSE()
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "totenergy" "-10.42800966 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "kinetic" "11.79557948 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "potential" "-22.22358913 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "eeenergy" "-1.30357077 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "ionion" "-12.77566756 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "localecp" "-9.17380267 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "nonlocalecp" "1.02945187 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "mpc" "-1.00378421 0.000001")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "flux" "-10.16222139 0.000001")
 ENDIF()

IF(NOT ENABLE_OFFLOAD)
  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj_kspace
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_kspace
                    det_qmc_short_kspace.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_KSPACE_SCALARS # VMC
                    )
ENDIF()

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.52336521 0.00001129")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "4.49057680 0.00001257")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-15.01394201 0.00000253")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.83205956 0.00000157")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77567016 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.00159178 0.00000147")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "0.93126037 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.40050428 0.00000225")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
 ELSE()
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.52335525 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "4.49059156 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-15.01394682 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.83205797 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77566756 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.00159889 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "0.93126167 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.40050042 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc
                    det_qmc_short_vmc_dmc.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_SCALARS # DMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-9.53268532 0.00001562")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "kinetic" "10.18368670 0.00001143")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "potential" "-19.71637202 0.00000463")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-2.38849840 0.00000125")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "ionion" "-12.77567016 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "localecp" "-4.96629415 0.0000047")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "0.41409069 0.000001 ")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "mpc" "-2.00359004 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "samples" "9 0.0")
 ELSE()
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-9.53268995 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "kinetic" "10.18368235 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "potential" "-19.71637230 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-2.38849743 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "localecp" "-4.96629755 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "0.41409012 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "mpc" "-2.00358802 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "samples" "9 0.0")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-dmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc_excited
                    det_qmc_short_vmc_dmc_excited.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_EXCITED_SCALARS # DMC
                    )


IF(QMC_MIXED_PRECISION)
  # DMC cycle 1
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "totenergy" "-9.53153130 0.0000139")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "kinetic" "9.00930034 0.00001647")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "potential" "-18.54083164 0.00000446")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "eeenergy" "-2.11983595 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "ionion" "-12.77567016 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "localecp" "-4.36659965 0.00000448")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "nonlocalecp" "0.72127412 0.00000178")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "mpc" "-1.93675985 0.00000122")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "samples" "9 0.0")
  # DMC cycle 2
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "totenergy" "-10.56917392 0.000027")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "kinetic" "14.41362921 0.00003329")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "potential" "-24.98280313 0.00003051")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "eeenergy" "-2.89879032 0.00000133")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "ionion" "-12.77567016 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "localecp" "-15.41387689 0.00001471")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "nonlocalecp" "6.10553425 0.00001448")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "mpc" "-2.61905734 0.00000179")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "samples" "9 0.0")
  # DMC cycle 3
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "totenergy" "-8.93215498 0.00001861")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "kinetic" "19.70827934 0.00003064")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "potential" "-28.64043432 0.00001851")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "eeenergy" "-1.15369571 0.00000776")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "ionion" "-12.77567016 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "localecp" "-16.18244873 0.00001536")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "nonlocalecp" "1.47138028 0.00000689")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "mpc" "-0.80671517 0.00000857")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "samples" "9 0.0")
 ELSE()
  # DMC cycle 1
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "totenergy" "-9.53151951 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "kinetic" "9.00930862 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "potential" "-18.54082813 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "eeenergy" "-2.11983454 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "localecp" "-4.36660098 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "nonlocalecp" "0.72127483 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "mpc" "-1.93675646 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "samples" "9 0.0")
  # DMC cycle 2
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "totenergy" "-10.56915192 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "kinetic" "14.41364745 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "potential" "-24.98279937 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "eeenergy" "-2.89879008 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "localecp" "-15.41389327 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "nonlocalecp" "6.10555142 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "mpc" "-2.61905530 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "samples" "9 0.0")
  # DMC cycle 3
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "totenergy" "-8.93215447 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "kinetic" "19.70828644 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "potential" "-28.64044091 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "eeenergy" "-1.15369726 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "localecp" "-16.18246211 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "nonlocalecp" "1.47138589 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "mpc" "-0.80671598 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "samples" "9 0.0")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-multi_dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_multi_dmc
                    det_qmc_short_vmc_multi_dmc.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_MULTI1_SCALARS # DMC cycle 1
                    2 DET_DIAMOND_DMC_MULTI2_SCALARS # DMC cycle 2
                    3 DET_DIAMOND_DMC_MULTI3_SCALARS # DMC cycle 3
                    )

IF(QMC_MIXED_PRECISION)
  # VMC cycle 1
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "totenergy" "-10.51638314 0.00000925")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "kinetic" "12.21076917 0.00000858")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "potential" "-22.72715231 0.00000257")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "eeenergy" "-1.39227245 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "ionion" "-12.77567016 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "localecp" "-8.27168446 0.0000035")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "nonlocalecp" "-0.28752524 0.00000126")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "mpc" "-0.95019080 0.00000112")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # VMC cycle 2
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "totenergy" "-10.56982090 0.00001796")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "kinetic" "13.81347725 0.0000143")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "potential" "-24.38329815 0.00000366")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "eeenergy" "-2.08653505 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "ionion" "-12.77567016  0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "localecp" "-10.59945535 0.00000897")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "nonlocalecp" "1.07836240 0.00000743")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "mpc" "-2.00098290 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "samples" "18 0.0")
  # DMC cycle 1 
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "totenergy" "-10.53594970 0.00000979")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "kinetic" "7.04049333 0.00000959")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "potential" "-17.57644302 0.00000443")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "eeenergy" "-1.15002939 0.00000344")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "ionion" "-12.77567016 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "localecp" "-3.78048184 0.00000184")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "nonlocalecp" "0.12973836 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "mpc" "-1.01499425 0.00000365")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # DMC cycle2
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "totenergy" "-9.72330854 0.0000083")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "kinetic" "11.45079701 0.00000847")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "potential" "-21.17410555 0.00000232")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "eeenergy" "-1.73937136 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "ionion" "-12.77567016 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "localecp" "-8.33546910 0.00000206")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "nonlocalecp" "1.67640507 0.00000168")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "mpc" "-1.50889162 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "samples" "18 0.0")
 ELSE()
  # VMC cycle 1
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "totenergy" "-10.51637540 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "kinetic" "12.21077547 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "potential" "-22.72715088 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "eeenergy" "-1.39227187 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "localecp" "-8.27169251 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "nonlocalecp" "-0.28751906 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "mpc" "-0.95018933 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # VMC cycle 2
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "totenergy" "-10.56981990 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "kinetic" "13.81348242 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "potential" "-24.38330232 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "eeenergy" "-2.08653334 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "localecp" "-10.59945467 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "nonlocalecp" "1.07835313 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "mpc" "-2.00098039 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "samples" "18 0.0")
  # DMC cycle 1
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "totenergy" "-10.53594951 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "kinetic" "7.04048475 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "potential" "-17.57643425 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "eeenergy" "-1.15002397 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "localecp" "-3.78048098 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "nonlocalecp" "0.12973814 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "mpc" "-1.01498843 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # DMC cycle 2
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "totenergy" "-9.72330413 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "kinetic" "11.45079979 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "potential" "-21.17410392 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "eeenergy" "-1.73936902 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "ionion" "-12.77566744 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "localecp" "-8.33547227 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "nonlocalecp" "1.67640481 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "mpc" "-1.50888780 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "samples" "18 0.0")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc-dmc-mwalkers_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc_mwalkers
                    det_qmc_short_vmc_dmc_mwalkers.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_VMC_MWALKERS1_SCALARS # VMC cycle 1
                    1 DET_DIAMOND_VMC_MWALKERS2_SCALARS # VMC cycle 2
                    2 DET_DIAMOND_DMC_MWALKERS1_SCALARS # DMC cycle 1
                    3 DET_DIAMOND_DMC_MWALKERS2_SCALARS # DMC cycle 2
                    )


IF((NOT QMC_MIXED_PRECISION) AND (NOT QMC_CUDA) AND (NOT QMC_COMPLEX))
  # DMC cycle 1
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "totenergy" "-11.15181221 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "kinetic" "11.34674025 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "potential" "-22.49855246 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "eeenergy" "-2.72967451 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "ionion" "-12.77566749 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "localecp" "-8.36105699 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "nonlocalecp" "1.36784653 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "mpc" "-2.66313760 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS1_SCALARS "samples" "16 0.0")
  # DMC cycle 2
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "totenergy" "-10.20986762 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "kinetic" "18.37238410 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "potential" "-28.58225171 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "eeenergy" "-3.11131238 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "ionion" "-12.77566749 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "localecp" "-11.73812555 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "nonlocalecp" "-0.95714628 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "mpc" "-3.02667676 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_LGTS2_SCALARS "samples" "18 0.0")

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc-dmc-multi_ts_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc_lg_ts
                    det_qmc_short_vmc_dmc_lg_ts.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_LGTS1_SCALARS # DMC cycle 1
                    2 DET_DIAMOND_DMC_LGTS2_SCALARS # DMC cycle 2
                    )
ENDIF()

  # traces tests
  set(traces_python_reqs numpy;h5py)
  CHECK_PYTHON_REQS(traces_python_reqs diamond-traces add_traces_tests)
  if (add_traces_tests)
   IF (QMC_MIXED_PRECISION)
    # traces with only scalars on
    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars.in.xml
      1 1
      check_traces.py -p qmc_trace_scalars -s '0,1' -m 'vmc,dmc' --dmc_steps_exclude=1 --tol=1e-6
      )
   ELSE()     
    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars.in.xml
      1 1
      check_traces.py -p qmc_trace_scalars -s '0,1' -m 'vmc,dmc' --dmc_steps_exclude=1
      )
   ENDIF()

   IF (QMC_MIXED_PRECISION)
    # traces with only selective scalars
    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars_sel
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars_selective.in.xml
      1 1
      check_traces.py -p qmc_trace_scalars_selective -s '0,1' -m 'vmc,dmc' -q 'Kinetic,ElecElec' --dmc_steps_exclude=1 --tol=1e-6
      )
   ELSE()
    SIMPLE_RUN_AND_CHECK(
      deterministic-diamondC_1x1x1_pp-vmc-dmc-trace_scalars_sel
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
      qmc_trace_scalars_selective.in.xml
      1 1
      check_traces.py -p qmc_trace_scalars_selective -s '0,1' -m 'vmc,dmc' -q 'Kinetic,ElecElec' --dmc_steps_exclude=1
      )
   ENDIF()
  endif()

# CUDA references
ELSE()
 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.55745923 0.0000213")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "9.01981971 0.00000356")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-19.57727588 0.00000255")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.98592540 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566732 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-5.21571393 0.00000221")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "0.40003076 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "mpc" "-1.72569333 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_SCALARS "flux" "-7.45176087 0.00000863")
 ELSE()
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.55746628 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "9.01980930 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-19.57727558 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.98592535 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566733 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-5.21571338 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "0.40003047 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "mpc" "-1.72569207 0.000001")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_SCALARS "flux" "-7.45173815 0.000001")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short
                    det_qmc_short.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_SCALARS # VMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "totenergy" "-9.95745569 0.00001281")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "kinetic" "10.21668907 0.00001385")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "potential" "-20.17414476 0.00000151")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "eeenergy" "-4.02104702 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "ionion" "-12.77566732 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "localecp" "-3.57226857 0.00000228")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "nonlocalecp" "0.19483814 0.00000166")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "mpc" "-3.62855170 0.0000011")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "flux" "-3.42024391 0.00002337")
 ELSE()
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "totenergy" "-9.95744357 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "kinetic" "10.21669957 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "potential" "-20.17414314 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "eeenergy" "-4.02104723 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "ionion" "-12.77566733 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "localecp" "-3.57227218 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "nonlocalecp" "0.19484358 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "mpc" "-3.62855318 0.000001")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_EXCITED_SCALARS "flux" "-3.42019506 0.000004")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_excited
                    det_qmc_short_excited.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_EXCITED_SCALARS # VMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.14372721 0.00006116")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "7.75142874 0.00005689")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-17.89515595 0.00000428")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-2.24186011 0.00000119")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77566732 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-2.90630818 0.00000161")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "0.02867966 0.00000147")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-1.96586493 0.00000228")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "9.51736092 0.00008676")
 ELSE()
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.14408211 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "7.75107010 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-17.89515221 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-2.24185969 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77566733 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-2.90630261 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "0.02867742 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-1.96586531 0.000001")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "9.51804679 0.000001")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short-meshf
                    det_qmc_short-meshf.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_MESHF_SCALARS # VMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.83933898 0.00000659")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "7.03523896 0.00000421")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-17.87457794 0.00000448")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.35584983 0.00000373")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77566732 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.50273251 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "-0.95202794 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.50931524 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
 ELSE()
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.83933516 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "7.03524341 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-17.87457857 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.35584811 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77566733 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.50273166 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "-0.95202769 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.50931484 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc
                    det_qmc_short_vmc_dmc.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_SCALARS # DMC
                    )

 IF(QMC_MIXED_PRECISION)
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-9.21905229 0.00000121")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "kinetic" "15.42550061 0.00000517")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "potential" "-24.64455290 0.00000412")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-3.60278140 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "ionion" "-12.77566732 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "localecp" "-5.26184234 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "-3.00426184 0.00000354")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "mpc" "-3.24225487 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "samples" "9 0.0")
 ELSE()
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "totenergy" "-9.21905527 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "kinetic" "15.42549492 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "potential" "-24.64455020 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "eeenergy" "-3.60278146 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "ionion" "-12.77566733 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "localecp" "-5.26184183 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "nonlocalecp" "-3.00425958 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "mpc" "-3.24225511 0.000001")
   LIST(APPEND DET_DIAMOND_DMC_EXCITED_SCALARS "samples" "9 0.0")
 ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-dmc_sdj_excited
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc_excited
                    det_qmc_short_vmc_dmc_excited.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_EXCITED_SCALARS # DMC
                    )

IF(QMC_MIXED_PRECISION)
  # DMC cycle 1
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "totenergy" "-10.07097496 0.0000068")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "kinetic" "16.06911449 0.00001353")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "potential" "-26.14008945 0.00000673")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "eeenergy" "-3.69351147 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "localecp" "-7.61944431 0.00000678")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "nonlocalecp" "-2.05146635 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "mpc" "-3.42701238 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "samples" "9 0.0")
  # DMC cycle 2
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "totenergy" "-10.09233269 0.00001359")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "kinetic" "13.20815324 0.00001392")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "potential" "-23.30048592 0.00000286")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "eeenergy" "-3.99940957 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "localecp" "-5.49404338 0.00000333")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "nonlocalecp" "-1.03136566 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "mpc" "-3.79251417 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "samples" "9 0.0")
  # DMC cycle 3
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "totenergy" "-10.50197534 0.0000106")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "kinetic" "12.25911588 0.0000084")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "potential" "-22.76109121 0.0000077")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "eeenergy" "-4.27297794 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "localecp" "-4.95181896 0.00000766")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "nonlocalecp" "-0.76062700 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "mpc" "-4.09025487 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "samples" "9 0.0")
ELSE()
  # DMC cycle 1
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "totenergy" "-10.07098874 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "kinetic" "16.06910421 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "potential" "-26.14009295 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "eeenergy" "-3.69351123 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "localecp" "-7.61944653 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "nonlocalecp" "-2.05146786 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "mpc" "-3.42701255 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI1_SCALARS "samples" "9 0.0")
  # DMC cycle 2
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "totenergy" "-10.09232889 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "kinetic" "13.20815760 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "potential" "-23.30048649 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "eeenergy" "-3.99940945 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "localecp" "-5.49404385 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "nonlocalecp" "-1.03136586 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "mpc" "-3.79251458 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI2_SCALARS "samples" "9 0.0")
  # DMC cycle 3
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "totenergy" "-10.50196932 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "kinetic" "12.25912048 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "potential" "-22.76108980 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "eeenergy" "-4.27297792 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "localecp" "-4.95181728 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "nonlocalecp" "-0.76062727 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "mpc" "-4.09025506 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MULTI3_SCALARS "samples" "9 0.0")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-multi_dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_multi_dmc
                    det_qmc_short_vmc_multi_dmc.in.xml
                    1 1
                    TRUE
                    1 DET_DIAMOND_DMC_MULTI1_SCALARS # DMC cycle 1
                    2 DET_DIAMOND_DMC_MULTI2_SCALARS # DMC cycle 2
                    3 DET_DIAMOND_DMC_MULTI3_SCALARS # DMC cycle 3
                    )

IF(QMC_MIXED_PRECISION)
  # VMC cycle 1
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "totenergy" "-11.11096324 0.00000291")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "kinetic" "11.87153330 0.00000486")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "potential" "-22.98249654 0.0000032")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "eeenergy" "-2.57496840 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "localecp" "-9.60549625 0.00000394")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "nonlocalecp" "1.97363542 0.00000156")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "mpc" "-2.33232140 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # VMC cycle 2
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "totenergy" "-10.85429855 0.00003797")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "kinetic" "11.77680157 0.00001599")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "potential" "-22.63110705 0.00000483")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "eeenergy" "-3.01702281 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "localecp" "-7.23401173 0.00000803")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "nonlocalecp" "0.39559601 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "mpc" "-2.80032609 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "samples" "18 0.0")
  # DMC cycle 1 
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "totenergy" "-10.31819204 0.00001753")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "kinetic" "15.62961326 0.00001111")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "potential" "-25.94780531 0.00000641")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "eeenergy" "-2.33834387 0.00000195")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "localecp" "-11.60746616 0.00000457")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "nonlocalecp" "0.77367224 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "mpc" "-2.03693208 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # DMC cycle 2 
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "totenergy" "-10.80471112 0.00000672")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "kinetic" "12.74023445 0.00001403")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "potential" "-23.54494494 0.00000795")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "eeenergy" "-3.38394953 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "ionion" "-12.77566732 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "localecp" "-8.70013857 0.0000099")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "nonlocalecp" "1.31481100 0.00000233")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "mpc" "-3.34110696 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "samples" "18 0.0")
 ELSE()
  # VMC cycle 1
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "totenergy" "-11.11096200 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "kinetic" "11.87153776 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "potential" "-22.98249977 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "eeenergy" "-2.57496810 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "localecp" "-9.60549471 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "nonlocalecp" "1.97363037 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "mpc" "-2.33232160 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # VMC cycle 2
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "totenergy" "-10.85431955 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "kinetic" "11.77679195 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "potential" "-22.63111150 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "eeenergy" "-3.01702362 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "localecp" "-7.23401682 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "nonlocalecp" "0.39559627 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "mpc" "-2.80032628 0.000001")
  LIST(APPEND DET_DIAMOND_VMC_MWALKERS2_SCALARS "samples" "18 0.0")
  # DMC cycle 1 
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "totenergy" "-10.31819817 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "kinetic" "15.62961311 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "potential" "-25.94781128 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "eeenergy" "-2.33834516 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "localecp" "-11.60747121 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "nonlocalecp" "0.77367241 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "mpc" "-2.03693292 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS1_SCALARS "samples" "18 0.0")
  # DMC cycle 2 
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "totenergy" "-10.80471716 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "kinetic" "12.74023623 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "potential" "-23.54495339 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "eeenergy" "-3.38394990 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "ionion" "-12.77566733 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "localecp" "-8.70014878 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "nonlocalecp" "1.31481262 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "mpc" "-3.34110760 0.000001")
  LIST(APPEND DET_DIAMOND_DMC_MWALKERS2_SCALARS "samples" "18 0.0")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc-dmc-mwalkers_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_vmc_dmc_mwalkers
                    det_qmc_short_vmc_dmc_mwalkers.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_VMC_MWALKERS1_SCALARS # VMC cycle 1
                    1 DET_DIAMOND_VMC_MWALKERS2_SCALARS # VMC cycle 2
                    2 DET_DIAMOND_DMC_MWALKERS1_SCALARS # DMC cycle 1
                    3 DET_DIAMOND_DMC_MWALKERS2_SCALARS # DMC cycle 2
                    )


ENDIF()
