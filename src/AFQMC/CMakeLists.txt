
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
  Hamiltonians/KPTHCHamiltonian.cpp
  Hamiltonians/HamiltonianFactory_Helper.cpp
  Hamiltonians/HamiltonianFactory.cpp
  Hamiltonians/THCHamiltonian.cpp
  )

IF(QMC_CUDA)
  SET(AFQMC_SRCS ${AFQMC_SRCS} 
    Memory/CUDA/cuda_utilities.cpp
    Kernels/determinant.cu 
    Kernels/adotpby.cu
    Kernels/fill_n.cu
    Kernels/uninitialized_fill_n.cu
    Kernels/uninitialized_copy_n.cu
    Kernels/axty.cu
    Kernels/adiagApy.cu
    Kernels/sum.cu
    Kernels/acAxpbB.cu
    Kernels/print.cu
    Kernels/setIdentity.cu
    Kernels/zero_complex_part.cu
    Kernels/batchedDot.cu
    Kernels/copy_n_cast.cu
    Kernels/ajw_to_waj.cu 
    Kernels/dot_wabn_wban.cu 
    Kernels/vKKwij_to_vwKiKj.cu 
    Kernels/KaKjw_to_QKajw.cu
    Kernels/vbias_from_v1.cu
    Kernels/KaKjw_to_KKwaj.cu
    Kernels/batched_dot_wabn_wban.cu
    Kernels/batched_Tab_to_Klr.cu
    Kernels/sampleGaussianRNG.cu
    Kernels/construct_X.cu
#    Kernels/myAtomicAdd.cu
  )
ENDIF(QMC_CUDA)

IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(afqmc ${AFQMC_SRCS})
  CUDA_ADD_CUBLAS_TO_TARGET(afqmc)
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
ENDIF()

