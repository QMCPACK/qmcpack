####################################
# QMC Drivers which implements various QMC algorithms
####################################
SET(QMCDRIVERS 
  SimpleFixedNodeBranch.cpp
  BranchIO.cpp
  QMCDriver.cpp
  QMCOptimize.cpp
  QMCChooseBestParameters.cpp
  QMCLinearOptimize.cpp
  QMCCorrelatedSamplingLinearOptimize.cpp
  QMCFixedSampleLinearOptimize.cpp
  QMCCSLinearOptimizeWFmanagerOMP.cpp
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
  ForwardWalking/FRSingleOMP.cpp
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
  ../Estimators/ReleasedNodeEnergyEstimator.cpp
  ../Estimators/AlternateReleasedNodeEnergyEstimator.cpp
  ../Estimators/WFMCOnlyEstimator.cpp
  ../Estimators/FWEstimator.cpp
  ../Estimators/LocalEnergyEstimatorHDF.cpp
  ../Estimators/EstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  ../Estimators/CollectablesEstimator.cpp
)

IF(NOT QMC_COMPLEX)
  SET(QMCDRIVERS  ${QMCDRIVERS}
    DMC/RNDMCOMP.cpp
    DMC/RNDMCUpdatePbyPCeperley.cpp
    DMC/RNDMCUpdatePbyPAlternate.cpp
    DMC/RNFactory.cpp
    )
ENDIF(NOT QMC_COMPLEX)


IF(QMC_CUDA)
   SET(QMCDRIVERS ${QMCDRIVERS}
   VMC/VMC_CUDA.cpp
   DMC/DMC_CUDA.cpp
   QMCCostFunctionCUDA.cpp)
ENDIF(QMC_CUDA)


IF(ENABLE_OPENMP) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
      VMC/VMCSingleOMP.cpp
      VMC/VMCLinearOptOMP.cpp
      WFMC/WFMCSingleOMP.cpp
      QMCCostFunctionOMP.cpp
     )
ENDIF(ENABLE_OPENMP) 

IF(HAVE_MPI) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 


IF(QMC_BUILD_LEVEL GREATER 1)
  SET(QMCDRIVERS ${QMCDRIVERS}
    PolymerEstimator.cpp
    MultiChain.cpp
    RQMCMultiple.cpp
    ../Estimators/ComboPolymerEstimator.cpp
    ../Estimators/CSEnergyEstimator.cpp
    ../Estimators/CSPolymerEstimator.cpp
    ../Estimators/MJPolymerEstimator.cpp
    ../Estimators/HFDHE2PolymerEstimator.cpp
  )
# REMOVE broken stuff
#    ../Estimators/RQMCMultipleEstimator.cpp
#    SpaceWarp.cpp
#    RQMCMultiWarp.cpp
#    VMC/VMCMultipleWarp.cpp
#    VMC/VMCPbyPMultiWarp.cpp
#    CorrelatedSampling/CSUpdateBase.cpp
#    CorrelatedSampling/CSVMC.cpp
#    CorrelatedSampling/CSVMCUpdateAll.cpp
#    CorrelatedSampling/CSVMCUpdatePbyP.cpp
#  IF(NOT QMC_COMPLEX)
#    SET(QMCDRIVERS ${QMCDRIVERS}
#    RQMCMultiplePbyP.cpp
#    )
#  ENDIF(NOT QMC_COMPLEX)

ENDIF(QMC_BUILD_LEVEL GREATER 1)


####################################
# create libqmc 
####################################
ADD_LIBRARY(qmcdriver ${QMCDRIVERS})
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver STATIC ${QMCDRIVERS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver SHARED ${QMCDRIVERS})
#ENDIF(QMC_BUILD_STATIC)
