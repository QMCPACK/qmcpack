
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

IF(ENABLE_SOA OR QMC_CUDA)
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_hybridrep_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_hybridrep.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping diamondC_1x1x1_hybridrep_pp-vmc_sdj test because neiher SOA nor CUDA is enabled")
ENDIF()

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short-meshf
                    qmc_short-meshf.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

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



# estimator tests
INCLUDE("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
set(stat_python_reqs numpy;h5py)
CHECK_PYTHON_REQS(stat_python_reqs diamond-estimator add_tests)
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
if (add_tests)
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

IF(QMC_MIXED_PRECISION)
  IF(QMC_COMPLEX)
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.65156025 0.00002215")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "16.03245759 0.00002646")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-26.68401784 0.00000661")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.64876057 0.00000199")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-13.45055710 0.00000658")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "1.19096740 0.00000788")
   LIST(APPEND DET_DIAMOND_SCALARS "mpc" "-1.45236706 0.0000041")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_SCALARS "flux" "-14.25663001 0.00009664")
  ELSE()
   LIST(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.65156025 0.00001842")
   LIST(APPEND DET_DIAMOND_SCALARS "kinetic" "16.03245759 0.00001793")
   LIST(APPEND DET_DIAMOND_SCALARS "potential" "-26.68401784 0.00000737")
   LIST(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.64876057 0.00000242")
   LIST(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_SCALARS "localecp" "-13.45055710 0.00001011")
   LIST(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "1.19096740 0.00000583")
   LIST(APPEND DET_DIAMOND_SCALARS "mpc" "-1.45236706 0.00000312")
   LIST(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_SCALARS "flux" "-14.25663001 0.00004515")
  ENDIF()
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
  IF(QMC_COMPLEX)
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.65617656 0.00024683")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "12.60316496 0.00025678")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-23.25934152 0.00000995")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-3.48893469 0.00000228")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-8.52555807 0.00001129")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "1.53081881 0.00000853")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-3.34182165 0.0000021")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "13.85780746 0.00068444")
  ELSE()
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "totenergy" "-10.65617656 0.00040281")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "kinetic" "12.60316496 0.00039482")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "potential" "-23.25934152 0.00001497")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "eeenergy" "-3.48893469 0.00000169")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "localecp" "-8.52555807 0.00002099")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "nonlocalecp" "1.53081881 0.00000969")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "mpc" "-3.34182165 0.00000308")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_MESHF_SCALARS "flux" "13.85780746 0.00089361")
  ENDIF()
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
  IF(QMC_COMPLEX)
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "totenergy" "-10.42800966 0.0000118")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "kinetic" "11.79557948 0.00001222")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "potential" "-22.22358913 0.00000605")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "eeenergy" "-1.30357077 0.00000482")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "localecp" "-9.17380267 0.00000821")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "nonlocalecp" "1.02945187 0.00000344")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "mpc" "-1.00378421 0.00000494")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "flux" "-10.16222139 0.00005138")
  ELSE()
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "totenergy" "-10.42800966 0.00000757")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "kinetic" "11.79557948 0.00000924")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "potential" "-22.22358913 0.000007")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "eeenergy" "-1.30357077 0.00000541")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "localecp" "-9.17380267 0.00000705")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "nonlocalecp" "1.02945187 0.00000206")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "mpc" "-1.00378421 0.00000653")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "samples" "9 0.0")
   LIST(APPEND DET_DIAMOND_KSPACE_SCALARS "flux" "-10.16222139 0.0001345")
  ENDIF()
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

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1_pp-vmc_sdj_kspace
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    det_qmc_short_kspace
                    det_qmc_short_kspace.in.xml
                    1 1
                    TRUE
                    0 DET_DIAMOND_KSPACE_SCALARS # VMC
                    )


IF(QMC_MIXED_PRECISION)
  IF(QMC_COMPLEX) 
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.52335525 0.00001161")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "4.49059156 0.00001792")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-15.01394682 0.0000074")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.83205797 0.0000032")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.00159889 0.0000081")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "0.93126167 0.00000158")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.40050042 0.0000054")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
  ELSE()
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "totenergy" "-10.52335525 0.0002301")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "kinetic" "4.49059156 0.00002106")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "potential" "-15.01394682 0.00000737")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "eeenergy" "0.83205797 0.0000032")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "ionion" "-12.77566756 0.00000265")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "localecp" "-4.00159889 0.00000809")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "nonlocalecp" "0.93126167 0.0000016")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "mpc" "1.40050042 0.0000054")
   LIST(APPEND DET_DIAMOND_DMC_SCALARS "samples" "9 0.0")
  ENDIF()
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
