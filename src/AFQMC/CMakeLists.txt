
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
#  Hamiltonians/KPTHCHamiltonian.cpp
  Hamiltonians/HamiltonianFactory_Helper.cpp
  Hamiltonians/HamiltonianFactory.cpp
  Hamiltonians/THCHamiltonian.cpp
  )

IF(QMC_CUDA)
  SET(AFQMC_SRCS ${AFQMC_SRCS} 
    Memory/CUDA/cuda_utilities.cpp
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
    Numerics/detail/CUDA/Kernels/sampleGaussianRNG.cu
    Numerics/detail/CUDA/Kernels/construct_X.cu
    Numerics/detail/CUDA/Kernels/reference_operations.cu
    Numerics/detail/CUDA/Kernels/term_by_term_matrix_vec.cu
#    Numerics/detail/CUDA/Kernels/myAtomicAdd.cu
  )
ENDIF(QMC_CUDA)

IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(afqmc ${AFQMC_SRCS})
  CUDA_ADD_CUBLAS_TO_TARGET(afqmc)
  target_link_libraries( afqmc ${CUDA_cusparse_LIBRARY} )
  target_link_libraries( afqmc ${CUDA_cusolver_LIBRARY} )
  target_link_libraries( afqmc ${CUDA_curand_LIBRARY} )
ELSE(QMC_CUDA)
  ADD_LIBRARY(afqmc ${AFQMC_SRCS})
ENDIF(QMC_CUDA)


IF (BUILD_UNIT_TESTS)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external_codes/catch)
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

