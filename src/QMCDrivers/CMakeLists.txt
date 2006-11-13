####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(DMCDRIVERS
  ../Estimators/LocalEnergyEstimator.cpp
  ../Estimators/ScalarEstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  WalkerControlBase.cpp
  DMC/DMCMoveAll.cpp
  DMC/DMCPbyP.cpp
  DMC/DMCFactory.cpp
  DMC/DMCUpdateBase.cpp
  DMC/DMCUpdateAll.cpp
  DMC/DMCUpdatePbyP.cpp
  DMC/DMCNonLocalUpdate.cpp
  DMC/WalkerControlFactory.cpp
  DMC/WalkerReconfiguration.cpp
)

IF(ENABLE_OPENMP) 
  SET(DMCDRIVERS ${DMCDRIVERS} 
  DMC/DMCPbyPOMP.cpp)
ENDIF(ENABLE_OPENMP) 

IF(HAVE_MPI) 
  SET(DMCDRIVERS ${DMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 


SET(QMCDRIVERS 
  SimpleFixedNodeBranch.cpp
  QMCDriver.cpp
  DummyQMC.cpp
  VMC.cpp
  VMCParticleByParticle.cpp
  ReptationMC.cpp
  QMCOptimize.cpp
  QMCCostFunction.cpp
  WaveFunctionTester.cpp
  ${DMCDRIVERS}
)

IF(NOT QMCPLUSPLUS_RELEASE)

  SET(QMCDRIVERS ${QMCDRIVERS}
    VMCMultiple.cpp
    VMCPbyPMultiple.cpp
    PolymerEstimator.cpp
    RQMCMultiple.cpp
    MultiChain.cpp
    SpaceWarp.cpp
    VMCMultipleWarp.cpp
    VMCPbyPMultiWarp.cpp
    RQMCMultiWarp.cpp
  )

#include GSLOptimize.cpp if libgsl is enabled
#  IF(HAVE_LIBGSL)
#    SET(QMCDRIVERS ${QMCDRIVERS} ../QMCDrivers/GSLOptimize.cpp)
#  ENDIF(HAVE_LIBGSL)
ENDIF(NOT QMCPLUSPLUS_RELEASE)

####################################
# create libqmc 
####################################
ADD_LIBRARY(qmcdriver ${QMCDRIVERS})
