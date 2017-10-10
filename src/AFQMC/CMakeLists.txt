
SET(AFQMC_DEBUG 3)

#-------------------------------------------------------------------
# Sources
#-------------------------------------------------------------------
SET (AFQMC_SRCS 
  AFQMCFactory.cpp
  Walkers/LocalWalkerHandler.cpp
  Walkers/DistWalkerHandler.cpp
  Drivers/AFQMCDriver.cpp 
  Drivers/selectedCI.cpp 
  Drivers/VMCDriver.cpp 
  Drivers/BenchmarkDriver.cpp
  Hamiltonians/SparseGeneralHamiltonian.cpp
  Wavefunctions/WavefunctionHandler.cpp
  Wavefunctions/PureSingleDeterminant.cpp
  Wavefunctions/MultiPureSingleDeterminant.cpp
  Wavefunctions/GeneralSingleDeterminant.cpp
  #Wavefunctions/MultiGeneralSingleDeterminant.cpp
  Propagators/phaseless_ImpSamp_ForceBias.cpp
  Propagators/VMCPropagator.cpp
  Numerics/SparseMatrixOperations.cpp
  Numerics/DenseMatrixOperations.cpp
  Sandbox/compare_libraries.cpp
  Wavefunctions/WavefunctionHelper.cpp
  Hamiltonians/DDProjector.cpp
  Hamiltonians/CCProjector.cpp
  Utilities/Utils.cpp
  ) 

ADD_LIBRARY(afqmc ${AFQMC_SRCS})

IF (BUILD_UNIT_TESTS)
  INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external_codes/catch)
  SUBDIRS(Matrix/tests)
  SUBDIRS(Numerics/tests)
ENDIF()

