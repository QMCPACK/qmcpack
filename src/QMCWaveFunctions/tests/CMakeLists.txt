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
    
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${QMCPACK_UNIT_TEST_DIR})

SET(SRC_DIR wavefunction)
SET(UTEST_DIR ${qmcpack_BINARY_DIR}/tests/${SRC_DIR})
#SET(UTEST_HDF_INPUT ${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT0 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT1 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT2 ${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae/pwscf.pwscf.h5)
SET(UTEST_HDF_INPUT3 ${qmcpack_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp/LiH-arb.pwscf.h5)
SET(UTEST_HDF_INPUT4 ${qmcpack_SOURCE_DIR}/tests/solids/monoO_noncollinear_1x1x1_pp/o2_45deg_spins.pwscf.h5)
SET(UTEST_HDF_INPUT5 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/C_diamond-twist-third.h5)
SET(UTEST_HDF_INPUT6 ${qmcpack_SOURCE_DIR}/src/QMCWaveFunctions/tests/lcao_spinor.h5)
SET(UTEST_HDF_INPUT7 ${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_MSD/LiH.orbs.h5)
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory "${UTEST_DIR}")
MAYBE_SYMLINK(${UTEST_HDF_INPUT0} ${UTEST_DIR}/diamondC_1x1x1.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT1} ${UTEST_DIR}/diamondC_2x1x1.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT2} ${UTEST_DIR}/bccH.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT3} ${UTEST_DIR}/LiH-arb.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT4} ${UTEST_DIR}/o2_45deg_spins.pwscf.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT5} ${UTEST_DIR}/C_diamond-twist-third.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT6} ${UTEST_DIR}/lcao_spinor.h5)
MAYBE_SYMLINK(${UTEST_HDF_INPUT7} ${UTEST_DIR}/LiH.orbs.h5)

SET(FILES_TO_COPY he_sto3g.wfj.xml ne_def2_svp.wfnoj.xml hcn.structure.xml hcn.wfnoj.xml hcn_downdet.cuspInfo.xml hcn_updet.cuspInfo.xml
    ethanol.structure.xml ethanol.wfnoj.xml ethanol_updet.cuspInfo.xml ethanol_downdet.cuspInfo.xml C_diamond-twist-third.structure.xml C_diamond-twist-third.wfj.xml cartesian_order.wfnoj.xml dirac_order.wfnoj.xml )

FOREACH(fname ${FILES_TO_COPY})
  EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${fname}" ${UTEST_DIR})
ENDFOREACH()

IF (NOT QMC_CUDA)
  IF (NOT QMC_COMPLEX)
    SET(MO_SRCS test_MO.cpp test_multiquintic_spline.cpp test_cartesian_ao.cpp)
    IF (NOT QMC_MIXED_PRECISION)
      SET(MO_SRCS ${MO_SRCS} test_soa_cusp_corr.cpp)
    ENDIF()
  ELSE()
      SET(MO_SRCS test_pyscf_complex_MO.cpp test_MO_spinor.cpp)
  ENDIF()
ENDIF()

SET(TRIALWF_SRC test_TrialWaveFunction.cpp test_wavefunction_factory.cpp
                test_TrialWaveFunction_diamondC_2x1x1.cpp test_TrialWaveFunction_He.cpp
                test_wavefunction_pool.cpp test_example_he.cpp test_lattice_gaussian.cpp)
SET(SPOSET_SRC test_spo_collection_input_spline.cpp test_spo_collection_input_LCAO_xml.cpp
               test_spo_collection_input_MSD_LCAO_h5.cpp test_einset.cpp
               test_einset_spinor.cpp test_hybridrep.cpp test_pw.cpp ${MO_SRCS})
SET(JASTROW_SRC test_bspline_jastrow.cpp test_counting_jastrow.cpp test_polynomial_eeI_jastrow.cpp
                test_rpa_jastrow.cpp test_user_jastrow.cpp test_kspace_jastrow.cpp test_pade_jastrow.cpp
                test_short_range_cusp_jastrow.cpp test_DiffTwoBodyJastrowOrbital.cpp)
SET(DETERMINANT_SRC FakeSPO.cpp makeRngSpdMatrix.cpp test_DiracDeterminantBatched.cpp test_DiracDeterminantBatched.cpp
                    test_multi_dirac_determinant.cpp test_dirac_matrix.cpp
                    test_multi_slater_determinant.cpp)

IF(ENABLE_CUDA)
  SET(DETERMINANT_SRC ${DETERMINANT_SRC} test_DiracMatrixComputeCUDA.cpp)
ENDIF(ENABLE_CUDA)
IF(ENABLE_OFFLOAD)
  SET(DETERMINANT_SRC ${DETERMINANT_SRC} test_DiracMatrixComputeOMPTarget.cpp)
ENDIF(ENABLE_OFFLOAD)



FOREACH(CATEGORY trialwf sposet jastrow determinant)
  SET(UTEST_EXE test_${SRC_DIR}_${CATEGORY})
  SET(UTEST_NAME deterministic-unit_${UTEST_EXE})
  STRING(TOUPPER "${CATEGORY}_SRC" SOURCE_FILE_VAR_NAME)
  ADD_EXECUTABLE(${UTEST_EXE} ${${SOURCE_FILE_VAR_NAME}})

  TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_main qmcwfs Math::BLAS_LAPACK Math::scalar_vector_functions utilities_for_test)
  IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} cuda_test_kernels qmcparticle qmcutil containers platform_omptarget)
  ENDIF()

  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
  SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})
  IF( "${ENABLE_SANITIZER}" STREQUAL "asan" ) 
    SET_TESTS_PROPERTIES(${UTEST_NAME} PROPERTIES LABELS noasan)
  ENDIF()
ENDFOREACH()

# do not merge with test_wavefunction_determinant until all of
# test_wavefunction_determinants tests can pass sanitizers
if (ENABLE_CUDA)
  SET(UTEST_EXE test_direct_inversion)
  SET(UTEST_NAME deterministic-unit_${UTEST_EXE})
  SET(DIRECT_INVERSION_SRC test_cuBLAS_LU.cpp)
  ADD_EXECUTABLE(${UTEST_EXE} test_cuBLAS_LU.cpp)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcwfs_direct_inversion_cuda utilities_for_test catch_main_no_mpi)
  IF(USE_OBJECT_TARGET)
    TARGET_LINK_LIBRARIES(${UTEST_EXE} platform_omptarget)
  ENDIF()
  TARGET_INCLUDE_DIRECTORIES(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
endif()

# benchmarking increases cost of test compilation so don't build all the tests with it.
if (ENABLE_CUDA)
  SET(UTEST_EXE benchmark_diracmatrixcompute)
  SET(UTEST_NAME deterministic-unit_${UTEST_EXE})
  SET(BENCHMARK_SRC makeRngSpdMatrix.cpp benchmark_DiracMatrixComputeCUDA.cpp)
  ADD_EXECUTABLE(${UTEST_EXE} ${BENCHMARK_SRC})
  TARGET_LINK_LIBRARIES(${UTEST_EXE} catch_benchmark_main qmcwfs Math::BLAS_LAPACK Math::scalar_vector_functions utilities_for_test)
  IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(${UTEST_EXE} qmcutil qmcparticle containers platform_omptarget)
  ENDIF()
  TARGET_INCLUDE_DIRECTORIES(${UTEST_EXE} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
  ADD_UNIT_TEST(${UTEST_NAME} 1 1 "${QMCPACK_UNIT_TEST_DIR}/${UTEST_EXE}")
endif()


