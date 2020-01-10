#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2019 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//                    Bryan Clark, bclark@Princeton.edu, Princeton University
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
  MCPopulation.cpp
  QMCDriver.cpp
  QMCDriverInput.cpp
  QMCDriverNew.cpp
  QMCOptimize.cpp
  QMCLinearOptimize.cpp
  QMCCorrelatedSamplingLinearOptimize.cpp
  QMCFixedSampleLinearOptimize.cpp
  Optimizers/DescentEngine.cpp
  Optimizers/HybridEngine.cpp
  QMCCostFunctionBase.cpp
  WaveFunctionTester.cpp
  WalkerControlBase.cpp
  CloneManager.cpp
  ContextForSteps.cpp
  Crowd.cpp
  QMCUpdateBase.cpp
  GreenFunctionModifiers/DriftModifierBuilder.cpp
  GreenFunctionModifiers/DriftModifierUNR.cpp
  VMC/VMCUpdatePbyP.cpp
  VMC/VMCUpdateAll.cpp
  VMC/VMCFactory.cpp
  VMC/VMCFactoryNew.cpp
  VMC/VMCDriverInput.cpp
  VMC/SOVMCUpdatePbyP.cpp
  DMC/DMC.cpp
  DMC/DMCUpdateAll.cpp
  DMC/DMCUpdatePbyPFast.cpp
  DMC/DMCFactory.cpp
  DMC/DMCFactoryNew.cpp
  DMC/DMCBatched.cpp
  DMC/DMCDriverInput.cpp
  DMC/WalkerControlFactory.cpp
  DMC/WalkerReconfiguration.cpp
  DMC/SODMCUpdatePbyPFast.cpp
  RMC/RMC.cpp
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
  ../Estimators/EstimatorManagerBase.cpp
  ../Estimators/EstimatorManagerCrowd.cpp
  ../Estimators/CollectablesEstimator.cpp
)

IF(QMC_CUDA)
   SET(QMCDRIVERS ${QMCDRIVERS}
   VMC/VMC_CUDA.cpp
   DMC/DMC_CUDA.cpp
   QMCCostFunctionCUDA.cpp)
ENDIF(QMC_CUDA)

SET(QMCDRIVERS ${QMCDRIVERS} 
    VMC/VMC.cpp
    VMC/VMCBatched.cpp
    VMC/VMCLinearOpt.cpp
    QMCCostFunction.cpp)

IF(HAVE_MPI) 
  SET(QMCDRIVERS ${QMCDRIVERS} 
  DMC/WalkerControlMPI.cpp
  DMC/WalkerReconfigurationMPI.cpp)
ENDIF(HAVE_MPI) 

# if using the LMYEngine interface, make sure we add the necessary functions
IF(BUILD_LMYENGINE_INTERFACE)
  SET(QMCDRIVERS ${QMCDRIVERS} 
      LMYEngineInterface/LMYE_QMCCostFunction.cpp
     )
ENDIF(BUILD_LMYENGINE_INTERFACE)

SET(QMCDRIVERS ${QMCDRIVERS}
  ../Estimators/CSEnergyEstimator.cpp
)


####################################
# create libqmc 
####################################
ADD_LIBRARY(qmcdriver ${QMCDRIVERS})
ADD_LIBRARY(qmcdriver_unit ${QMCDRIVERS})
USE_FAKE_RNG(qmcdriver_unit)
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver STATIC ${QMCDRIVERS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcdriver SHARED ${QMCDRIVERS})
#ENDIF(QMC_BUILD_STATIC)
