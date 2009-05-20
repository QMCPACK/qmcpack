####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(QMCDRIVERS 
  SimpleFixedNodeBranch.cpp
  BranchIO.cpp
  QMCDriver.cpp
  QMCOptimize.cpp
  QMCLinearOptimize.cpp
  QMCCostFunctionBase.cpp
  QMCCostFunctionSingle.cpp
  WaveFunctionTester.cpp
  WalkerControlBase.cpp
  CloneManager.cpp
  QMCUpdateBase.cpp
  VMC/VMCUpdatePbyP.cpp
  VMC/VMCUpdateAll.cpp
  VMC/VMCSingle.cpp
  ForwardWalking/FWSingleOMP.cpp
  ForwardWalking/FWSingleMPI.cpp
  VMC/VMCFactory.cpp
  WFMC/WFMCUpdateAll.cpp
  DMC/DMCOMP.cpp
  DMC/DMCUpdateAll.cpp
  DMC/DMCUpdatePbyP.cpp
  DMC/DMCUpdatePbyPFast.cpp
  DMC/DMCFactory.cpp
  DMC/WalkerControlFactory.cpp
  DMC/WalkerReconfiguration.cpp
  ../Estimators/LocalEnergyEstimator.cpp
  ../Estimators/WFMCOnlyEstimator.cpp
  ../Estimators/FWEstimator.cpp
  ../Estimators/LocalEnergyEstimatorHDF.cpp
  ../Estimators/EstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  ../Estimators/CollectablesEstimator.cpp
)

IF(ENABLE_OPENMP) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
      VMC/VMCSingleOMP.cpp
      WFMC/WFMCSingleOMP.cpp
      QMCCostFunctionOMP.cpp
     )
ENDIF(ENABLE_OPENMP) 

IF(HAVE_MPI) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 


IF(QMC_BUILD_COMPLETE)
  SET(QMCDRIVERS ${QMCDRIVERS}
    PolymerEstimator.cpp
    RQMCMultiple.cpp
    RQMCMultiplePbyP.cpp
    MultiChain.cpp
    SpaceWarp.cpp
    RQMCMultiWarp.cpp
    VMC/VMCMultipleWarp.cpp
    VMC/VMCPbyPMultiWarp.cpp
    ../Estimators/RQMCMultipleEstimator.cpp
    ../Estimators/ComboPolymerEstimator.cpp
    ../Estimators/CSEnergyEstimator.cpp
    ../Estimators/CSPolymerEstimator.cpp
    ../Estimators/MJPolymerEstimator.cpp
    ../Estimators/HFDHE2PolymerEstimator.cpp
    CorrelatedSampling/CSUpdateBase.cpp
    CorrelatedSampling/CSVMC.cpp
    CorrelatedSampling/CSVMCUpdateAll.cpp
    CorrelatedSampling/CSVMCUpdatePbyP.cpp
  )
ENDIF(QMC_BUILD_COMPLETE)

####################################
# create libqmc 
####################################
ADD_LIBRARY(qmcdriver ${QMCDRIVERS})
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver STATIC ${QMCDRIVERS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver SHARED ${QMCDRIVERS})
#ENDIF(QMC_BUILD_STATIC)
