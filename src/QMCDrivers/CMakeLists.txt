####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(DMCDRIVERS
  ../Estimators/LocalEnergyEstimator.cpp
  ../Estimators/ScalarEstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  ../Estimators/RQMCMultipleEstimator.cpp
  WalkerControlBase.cpp
  CloneManager.cpp
  QMCUpdateBase.cpp
  VMC/VMCUpdatePbyP.cpp
  VMC/VMCUpdateAll.cpp
  VMC/VMCSingle.cpp
  VMC/VMCSingleOMP.cpp
  VMC/VMCFactory.cpp
  DMC/DMCUpdateAll.cpp
  DMC/DMCUpdatePbyP.cpp
  DMC/DMCNonLocalUpdate.cpp
  DMC/DMC.cpp
  DMC/DMCFactory.cpp
  DMC/WalkerControlFactory.cpp
  DMC/WalkerReconfiguration.cpp
)

IF(ENABLE_OPENMP) 
  SET(DMCDRIVERS ${DMCDRIVERS} DMC/DMCOMP.cpp)
ENDIF(ENABLE_OPENMP) 

IF(HAVE_MPI) 
  SET(DMCDRIVERS ${DMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 


SET(QMCDRIVERS 
  SimpleFixedNodeBranch.cpp
  QMCDriver.cpp
  QMCOptimize.cpp
  QMCCostFunction.cpp
  WaveFunctionTester.cpp
  ${DMCDRIVERS}
)

IF(NOT QMCPLUSPLUS_RELEASE)

  SET(QMCDRIVERS ${QMCDRIVERS}
    PolymerEstimator.cpp
    RQMCMultiple.cpp
    MultiChain.cpp
    SpaceWarp.cpp
    RQMCMultiWarp.cpp
    VMC/VMCMultiple.cpp
    VMC/VMCPbyPMultiple.cpp
    VMC/VMCMultipleWarp.cpp
    VMC/VMCPbyPMultiWarp.cpp
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
