#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
#//////////////////////////////////////////////////////////////////////////////////////
    
    


####################################
# QMC Drivers which implements various QMC algorithms
####################################

SET(QMCDRIVERS 
  SimpleFixedNodeBranch.cpp
  BranchIO.cpp
  QMCDriver.cpp
  QMCOptimize.cpp
  QMCLinearOptimize.cpp
  QMCCorrelatedSamplingLinearOptimize.cpp
  QMCFixedSampleLinearOptimize.cpp
  QMCCSLinearOptimizeWFmanagerOMP.cpp
  QMCCostFunctionBase.cpp
  WaveFunctionTester.cpp
  WalkerControlBase.cpp
  CloneManager.cpp
  QMCUpdateBase.cpp
  VMC/VMCUpdatePbyP.cpp
  VMC/VMCUpdateAll.cpp
  VMC/VMCFactory.cpp
  DMC/DMCOMP.cpp
  DMC/DMCUpdateAll.cpp
  DMC/DMCUpdatePbyPFast.cpp
  DMC/DMCFactory.cpp
  DMC/WalkerControlFactory.cpp
  DMC/WalkerReconfiguration.cpp
  RMC/RMCSingleOMP.cpp
  RMC/RMCUpdatePbyP.cpp
  RMC/RMCUpdateAll.cpp
  RMC/RMCFactory.cpp
  CorrelatedSampling/CSVMC.cpp
  CorrelatedSampling/CSVMCUpdateAll.cpp
  CorrelatedSampling/CSVMCUpdatePbyP.cpp
  CorrelatedSampling/CSUpdateBase.cpp
  ../Estimators/CSEnergyEstimator.cpp
  ../Estimators/LocalEnergyEstimator.cpp
  ../Estimators/RMCLocalEnergyEstimator.cpp
  ../Estimators/LocalEnergyEstimatorHDF.cpp
  ../Estimators/EstimatorManagerBase.cpp
#  ../Estimators/EstimatorManager.cpp
  ../Estimators/MultipleEnergyEstimator.cpp
  ../Estimators/CollectablesEstimator.cpp
)

IF(QMC_CUDA)
   SET(QMCDRIVERS ${QMCDRIVERS}
   VMC/VMC_CUDA.cpp
   DMC/DMC_CUDA.cpp
   QMCCostFunctionCUDA.cpp)
ENDIF(QMC_CUDA)

  SET(QMCDRIVERS ${QMCDRIVERS} 
      VMC/VMCSingleOMP.cpp
      VMC/VMCLinearOptOMP.cpp
      QMCCostFunctionOMP.cpp
     )

IF(HAVE_MPI) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 

# if using the LMYEngine interface, make sure we add the necessary functions
IF(BUILD_LMYENGINE_INTERFACE)
  SET(QMCDRIVERS ${QMCDRIVERS} 
      LMYEngineInterface/LMYE_QMCCostFunctionOMP.cpp
     )
ENDIF(BUILD_LMYENGINE_INTERFACE)

IF(QMC_BUILD_LEVEL GREATER 1)
  SET(QMCDRIVERS ${QMCDRIVERS}
#     PolymerEstimator.cpp
#     MultiChain.cpp
#     RQMCMultiple.cpp
    ../Estimators/CSEnergyEstimator.cpp
#     ../Estimators/CSPolymerEstimator.cpp
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
IF (HAS_TARGET_COMPILE_DEFINITIONS)
  ADD_LIBRARY(qmcdriver_unit ${QMCDRIVERS})
  USE_FAKE_RNG(qmcdriver_unit)
ENDIF()
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver STATIC ${QMCDRIVERS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver SHARED ${QMCDRIVERS})
#ENDIF(QMC_BUILD_STATIC)
