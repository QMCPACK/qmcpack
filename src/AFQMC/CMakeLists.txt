SET(AFQMC_DEBUG 3)

#-------------------------------------------------------------------
# Sources
#-------------------------------------------------------------------
SET (AFQMC_SRCS
  AFQMCFactory.cpp
  Drivers/DriverFactory.cpp
  Drivers/AFQMCDriver.cpp
#  Drivers/BenchmarkDriver.cpp
  Propagators/AFQMCBasePropagator.cpp
  Propagators/PropagatorFactory.cpp
  Wavefunctions/WavefunctionFactory.cpp
  Utilities/readWfn.cpp
  Hamiltonians/HSPotential_Helpers.cpp
  Hamiltonians/FactorizedSparseHamiltonian.cpp
  Hamiltonians/KPFactorizedHamiltonian.cpp
  Hamiltonians/RealDenseHamiltonian.cpp
  Hamiltonians/RealDenseHamiltonian_v2.cpp
#  Hamiltonians/KPTHCHamiltonian.cpp
  Hamiltonians/HamiltonianFactory_Helper.cpp
  Hamiltonians/HamiltonianFactory.cpp
  Hamiltonians/THCHamiltonian.cpp
  )

IF(ENABLE_CUDA)
  SET(AFQMC_SRCS ${AFQMC_SRCS} 
    Memory/CUDA/cuda_utilities.cpp
    Memory/CUDA/cuda_init.cpp
    Numerics/detail/CUDA/Kernels/determinant.cu 
    Numerics/detail/CUDA/Kernels/adotpby.cu
    Numerics/detail/CUDA/Kernels/fill_n.cu
    Numerics/detail/CUDA/Kernels/uninitialized_fill_n.cu
    Numerics/detail/CUDA/Kernels/uninitialized_copy_n.cu
    Numerics/detail/CUDA/Kernels/axty.cu
    Numerics/detail/CUDA/Kernels/adiagApy.cu
    Numerics/detail/CUDA/Kernels/sum.cu
    Numerics/detail/CUDA/Kernels/acAxpbB.cu
    Numerics/detail/CUDA/Kernels/print.cu
    Numerics/detail/CUDA/Kernels/setIdentity.cu
    Numerics/detail/CUDA/Kernels/zero_complex_part.cu
    Numerics/detail/CUDA/Kernels/batchedDot.cu
    Numerics/detail/CUDA/Kernels/copy_n_cast.cu
    Numerics/detail/CUDA/Kernels/ajw_to_waj.cu 
    Numerics/detail/CUDA/Kernels/vKKwij_to_vwKiKj.cu 
    Numerics/detail/CUDA/Kernels/KaKjw_to_QKajw.cu
    Numerics/detail/CUDA/Kernels/vbias_from_v1.cu
    Numerics/detail/CUDA/Kernels/KaKjw_to_KKwaj.cu
    Numerics/detail/CUDA/Kernels/batched_dot_wabn_wban.cu
    Numerics/detail/CUDA/Kernels/batched_Tab_to_Klr.cu
    Numerics/detail/CUDA/Kernels/dot_wabn.cu
    Numerics/detail/CUDA/Kernels/Tab_to_Kl.cu
    Numerics/detail/CUDA/Kernels/sampleGaussianRNG.cu
    Numerics/detail/CUDA/Kernels/construct_X.cu
    Numerics/detail/CUDA/Kernels/reference_operations.cu
    Numerics/detail/CUDA/Kernels/term_by_term_matrix_vec.cu
    Numerics/detail/CUDA/Kernels/axpyBatched.cu
#    Numerics/detail/CUDA/Kernels/myAtomicAdd.cu
  )
ENDIF(ENABLE_CUDA)

IF(ENABLE_CUDA)
  CUDA_ADD_LIBRARY(afqmc ${AFQMC_SRCS})
  CUDA_ADD_CUBLAS_TO_TARGET(afqmc)
  TARGET_LINK_LIBRARIES( afqmc PRIVATE ${CUDA_cusparse_LIBRARY} )
  TARGET_LINK_LIBRARIES( afqmc PRIVATE ${CUDA_cusolver_LIBRARY} )
  TARGET_LINK_LIBRARIES( afqmc PRIVATE ${CUDA_curand_LIBRARY} )
ELSE(ENABLE_CUDA)
  ADD_LIBRARY(afqmc ${AFQMC_SRCS})
ENDIF(ENABLE_CUDA)

TARGET_LINK_LIBRARIES(afqmc PUBLIC qmcutil)
TARGET_LINK_LIBRARIES(afqmc PRIVATE Math::BLAS_LAPACK Boost::boost)
IF(HAVE_MPI)
  TARGET_COMPILE_DEFINITIONS(afqmc PRIVATE "USE_MPI")
ENDIF(HAVE_MPI)

IF (BUILD_UNIT_TESTS)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external_codes/catch)
  SET(AFQMC_UNIT_TEST_INPUT_NOSYM
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/ham_chol_sc.h5
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/wfn_rhf.dat
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/ham_chol_sc.h5
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/wfn_uhf.dat
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/ham_chol_sc.h5
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/wfn_msd.dat
      #${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/ham_chol_sc.h5
      #${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_dzvp/wfn_phmsd.dat
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_szv/ham_thc_sc.h5
      ${qmcpack_SOURCE_DIR}/tests/afqmc/C_1x1x1_szv/wfn_rhf.dat
      )
  SET(AFQMC_UNIT_TEST_INPUT_SYM
      ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/ham_chol_uc.h5
      ${qmcpack_SOURCE_DIR}/tests/afqmc/He_2x2x2_dzv/wfn_rhf.dat
      )
  SET(AFQMC_UNIT_TEST_INPUT_REAL
      ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/ham_chol.h5
      ${qmcpack_SOURCE_DIR}/tests/afqmc/Ne_cc-pvdz/wfn_rhf.dat
      )
  IF (QMC_COMPLEX)
    IF (ENABLE_CUDA)
      SET(AFQMC_UNIT_TEST_INPUTS ${AFQMC_UNIT_TEST_INPUT_SYM})
    ELSE()
      SET(AFQMC_UNIT_TEST_INPUTS
          ${AFQMC_UNIT_TEST_INPUT_SYM}
          ${AFQMC_UNIT_TEST_INPUT_NOSYM}
          )
    ENDIF()
  ELSE()
    SET(AFQMC_UNIT_TEST_INPUTS ${AFQMC_UNIT_TEST_INPUT_REAL})
  ENDIF()
  SUBDIRS(Matrix/tests)
  SUBDIRS(Numerics/tests)
  #SUBDIRS(Numerics/shm_tests)
  SUBDIRS(SlaterDeterminantOperations/tests)
  SUBDIRS(Walkers/tests)
  SUBDIRS(Hamiltonians/tests)
  SUBDIRS(HamiltonianOperations/tests)
  SUBDIRS(Wavefunctions/tests)
  SUBDIRS(Propagators/tests)
  SUBDIRS(Estimators/tests)
ENDIF()

