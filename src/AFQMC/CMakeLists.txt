
SET(AFQMC_DEBUG 3)

#-------------------------------------------------------------------
# Sources
#-------------------------------------------------------------------
SET (AFQMC_SRCS
  AFQMCFactory.cpp
  Walkers/SharedWalkerSet.cpp
  Drivers/DriverFactory.cpp
  Drivers/AFQMCDriver.cpp
#  Drivers/BenchmarkDriver.cpp
  Propagators/AFQMCSharedPropagator.cpp
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

ADD_LIBRARY(afqmc ${AFQMC_SRCS})

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

