
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.491422243970 0.001977730412")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.434445621324 0.016479005831")
  LIST(APPEND DIAMOND_SCALARS "potential" "-21.925867865291 0.016428022009")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.700464664942 0.004217783716")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.775667433000 0.001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.047030147668 0.018491182454")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.597294379836 0.005597158667")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
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

  LIST(APPEND DIAMOND_KSPACE_SCALARS "totenergy" "-10.500668187915 0.001729410454")
  LIST(APPEND DIAMOND_KSPACE_SCALARS "variance"  "0.315109978771 0.064518513918")
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
  LIST(APPEND DIAMOND_DMC_SCALARS "totenergy"   "-10.531304509572 0.002448453821")
  LIST(APPEND DIAMOND_DMC_SCALARS "kinetic"      "11.484139015626 0.027421109472")
  LIST(APPEND DIAMOND_DMC_SCALARS "potential"   "-22.015443525201 0.027947337747")
  LIST(APPEND DIAMOND_DMC_SCALARS "localecp"     "-7.150962399967 0.035541672643")
  LIST(APPEND DIAMOND_DMC_SCALARS "nonlocalecp"  "0.627729831525 0.008032960873")
  LIST(APPEND DIAMOND_DMC_SCALARS "eeenergy"    "-2.716543566587 0.007869453952")

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

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.491438986865 0.000065863965")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "122880000 0.0")
  LIST(APPEND LONG_DIAMOND_SCALARS "flux" "0.0 0.03")

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
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "totenergy" "-10.491484837376 0.003952892737")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "kinetic" "11.434327164597 0.038998414966")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "potential" "-21.925812001976 0.039432901659")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "eeenergy" "-2.700608037894 0.011063923511")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "ionion" "-12.775667436000 0.0004")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "localecp" "-7.046553570494 0.048384004511")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "nonlocalecp" "0.597017042843 0.011269004144")
  LIST(APPEND LONG_DIAMOND_ALLP_SCALARS "samples" "1280000 0.0")
  #DMC ref
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "totenergy"   "-10.531001135281 0.000809439354")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "kinetic"      "11.486688574567 0.008932748191")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "potential"   "-22.017689697854 0.008806258877")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "localecp"     "-7.152139490279 0.010839755811")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "nonlocalecp"  "0.627051485859 0.002550989598")
  LIST(APPEND LONG_DIAMOND_DMC_ALLP_SCALARS "eeenergy"    "-2.716934257006 0.002460783830")

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

  LIST(APPEND LONG_DIAMOND_KSPACE_SCALARS "totenergy" "-10.500668187915 0.000547271940")
  LIST(APPEND LONG_DIAMOND_KSPACE_SCALARS "variance"  "0.315109978771 0.020416883784")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj_kspace
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_kspace
                    qmc_long_kspace.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_KSPACE_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.531304509572 0.000775990998")
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
