#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


set(SRC_DIR wavefunction)
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})
#SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp/pwscf.pwscf.h5)
set(UTEST_HDF_INPUT0 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
set(UTEST_HDF_INPUT1 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp/pwscf.pwscf.h5)
set(UTEST_HDF_INPUT2 ${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae/pwscf.pwscf.h5)
set(UTEST_HDF_INPUT3 ${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/LiH-arb.pwscf.h5)
set(UTEST_HDF_INPUT4 ${qmcpack_SOURCE_DIR}/tests/solids/monoO_noncollinear_1x1x1_pp/o2_45deg_spins.pwscf.h5)
set(UTEST_HDF_INPUT5 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/C_diamond-twist-third.h5)
set(UTEST_HDF_INPUT6 ${qmcpack_SOURCE_DIR}/src/QMCWaveFunctions/tests/lcao_spinor.h5)
set(UTEST_HDF_INPUT7 ${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD/LiH.orbs.h5)
set(UTEST_HDF_INPUT8 ${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD/LiH.Multidet.h5)
set(UTEST_HDF_INPUT9 ${qmcpack_SOURCE_DIR}/tests/converter/test_Bi_dirac/gold.orbs.h5)

maybe_symlink(${UTEST_HDF_INPUT0} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT1} ${UTEST_DIR}/diamondC_2x1x1.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT2} ${UTEST_DIR}/bccH.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT3} ${UTEST_DIR}/LiH-arb.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT4} ${UTEST_DIR}/o2_45deg_spins.pwscf.h5)
maybe_symlink(${UTEST_HDF_INPUT5} ${UTEST_DIR}/C_diamond-twist-third.h5)
maybe_symlink(${UTEST_HDF_INPUT6} ${UTEST_DIR}/lcao_spinor.h5)
maybe_symlink(${UTEST_HDF_INPUT7} ${UTEST_DIR}/LiH.orbs.h5)
maybe_symlink(${UTEST_HDF_INPUT8} ${UTEST_DIR}/LiH.Multidet.h5)
maybe_symlink(${UTEST_HDF_INPUT9} ${UTEST_DIR}/Bi.orbs.h5)

set(FILES_TO_COPY
    he_sto3g.wfj.xml
    ne_def2_svp.wfnoj.xml
    hcn.structure.xml
    hcn.wfnoj.xml
    hcn_downdet.cuspInfo.xml
    hcn_updet.cuspInfo.xml
    ethanol.structure.xml
    ethanol.wfnoj.xml
    ethanol_updet.cuspInfo.xml
    ethanol_downdet.cuspInfo.xml
    C_diamond-twist-third.structure.xml
    C_diamond-twist-third.wfj.xml
    cartesian_order.wfnoj.xml
    dirac_order.wfnoj.xml)

foreach(fname ${FILES_TO_COPY})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${fname}" ${UTEST_DIR})
endforeach()

if(NOT QMC_CUDA)
  if(NOT QMC_COMPLEX)
    set(MO_SRCS test_MO.cpp test_multiquintic_spline.cpp test_cartesian_ao.cpp)
    if(NOT QMC_MIXED_PRECISION)
      set(MO_SRCS ${MO_SRCS} test_soa_cusp_corr.cpp)
    endif()
  else()
    set(MO_SRCS test_pyscf_complex_MO.cpp test_MO_spinor.cpp)
  endif()
endif()

set(COMMON_SRC test_variable_set.cpp)
set(TRIALWF_SRC
    test_TrialWaveFunction.cpp
    test_wavefunction_factory.cpp
    test_TrialWaveFunction_diamondC_2x1x1.cpp
    test_TrialWaveFunction_He.cpp
    test_wavefunction_pool.cpp
    test_example_he.cpp
    test_lattice_gaussian.cpp
    test_TWFGrads.cpp)
set(SPOSET_SRC
    test_spo_collection_input_spline.cpp
    test_spo_collection_input_LCAO_xml.cpp
    test_spo_collection_input_MSD_LCAO_h5.cpp
    test_einset.cpp
    test_einset_spinor.cpp
    test_CompositeSPOSet.cpp
    test_hybridrep.cpp
    test_pw.cpp
    ${MO_SRCS})
set(JASTROW_SRC
    test_bspline_jastrow.cpp
    test_counting_jastrow.cpp
    test_polynomial_eeI_jastrow.cpp
    test_rpa_jastrow.cpp
    test_user_jastrow.cpp
    test_kspace_jastrow.cpp
    test_pade_jastrow.cpp
    test_short_range_cusp_jastrow.cpp
    test_J1OrbitalSoA.cpp
    test_J1Spin.cpp
    test_J2_bspline.cpp
    test_J2_derivatives.cpp)
set(DETERMINANT_SRC
    FakeSPO.cpp
    test_DiracDeterminant.cpp
    test_DiracDeterminantBatched.cpp
    test_multi_dirac_determinant.cpp
    test_dirac_matrix.cpp
    test_ci_configuration.cpp
    test_multi_slater_determinant.cpp)

if(ENABLE_CUDA)
  set(DETERMINANT_SRC ${DETERMINANT_SRC} test_DiracMatrixComputeCUDA.cpp)
endif(ENABLE_CUDA)
if(ENABLE_OFFLOAD)
  set(DETERMINANT_SRC ${DETERMINANT_SRC} test_DiracMatrixComputeOMPTarget.cpp)
endif(ENABLE_OFFLOAD)

foreach(CATEGORY common trialwf sposet jastrow determinant)
  set(UTEST_EXE test_${SRC_DIR}_${CATEGORY})
  set(UTEST_NAME deterministic-unit_${UTEST_EXE})
  string(TOUPPER "${CATEGORY}_SRC" SOURCE_FILE_VAR_NAME)
  add_executable(${UTEST_EXE} ${${SOURCE_FILE_VAR_NAME}})

  target_link_libraries(${UTEST_EXE} catch_main qmcwfs platform_LA platform_runtime
                        utilities_for_test container_testing)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} qmcparticle qmcutil platform_omptarget_LA)
  endif()

  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
  set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
endforeach()

# do not merge with test_wavefunction_determinant until all of
# test_wavefunction_determinants tests can pass sanitizers
if(ENABLE_CUDA)
  set(UTEST_EXE test_direct_inversion)
  set(UTEST_NAME deterministic-unit_${UTEST_EXE})
  set(DIRECT_INVERSION_SRC test_cuBLAS_LU.cpp)
  add_executable(${UTEST_EXE} test_cuBLAS_LU.cpp)
  target_link_libraries(${UTEST_EXE} qmcwfs_direct_inversion_cuda utilities_for_test catch_main_no_mpi)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} platform_omptarget_LA)
  endif()
  target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
endif()

if(ENABLE_CUDA AND BUILD_MICRO_BENCHMARKS)
  set(UTEST_EXE benchmark_diracmatrixcompute)
  set(UTEST_NAME deterministic-unit_${UTEST_EXE})
  set(BENCHMARK_SRC benchmark_DiracMatrixComputeCUDA.cpp)
  add_executable(${UTEST_EXE} ${BENCHMARK_SRC})
  target_link_libraries(${UTEST_EXE} catch_main qmcwfs platform_LA platform_runtime
                        utilities_for_test container_testing)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${UTEST_EXE} qmcutil qmcparticle platform_omptarget_LA)
  endif()
#  target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../..")
  target_include_directories(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
  add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
endif()
