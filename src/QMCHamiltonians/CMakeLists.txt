#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
#//                    D.C. Yang, University of Illinois at Urbana-Champaign
#//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////
    
    
SET(HAMSRCS
  HamiltonianPool.cpp
  ACForce.cpp
  OperatorBase.cpp
  QMCHamiltonian.cpp
  Hdispatcher.cpp
  CoulombPBCAA.cpp
  CoulombPBCAB.cpp
  NonLocalTOperator.cpp
  ForwardWalking.cpp
  PairCorrEstimator.cpp
  DensityEstimator.cpp
  SkPot.cpp
  SkEstimator.cpp
  SkAllEstimator.cpp
  MomentumEstimator.cpp
  ForceBase.cpp
  HamiltonianFactory.cpp
  CoulombPotentialFactory.cpp
  HarmonicExternalPotential.cpp
  GridExternalPotential.cpp
  StaticStructureFactor.cpp
  OrbitalImages.cpp
  SpinDensity.cpp
  SpeciesKineticEnergy.cpp
  LatticeDeviationEstimator.cpp
  SpaceWarpTransformation.cpp
  )

IF(OHMMS_DIM MATCHES 3)
  SET(HAMSRCS ${HAMSRCS}
    ChiesaCorrection.cpp
    ForceCeperley.cpp
    ForceChiesaPBCAA.cpp
    StressPBC.cpp
    LocalECPotential.cpp
    NonLocalECPComponent.cpp
    NonLocalECPotential.deriv.cpp
    NonLocalECPotential.cpp
    L2Potential.cpp
    SOECPComponent.cpp
    SOECPotential.cpp
    ECPotentialBuilder.cpp
    ECPComponentBuilder.cpp
    ECPComponentBuilder.1.cpp
    ECPComponentBuilder.2.cpp
    ECPComponentBuilder_L2.cpp
    )


  IF(NOT REMOVE_TRACEMANAGER)
    SET(HAMSRCS ${HAMSRCS}
      ../Estimators/TraceManager.cpp
      ReferencePoints.cpp
      SpaceGrid.cpp
      EnergyDensityEstimator.cpp
      DensityMatrices1B.cpp
      )
  ENDIF()

  IF(HAVE_LIBFFTW)
    SET(HAMSRCS ${HAMSRCS} MPC.cpp)
  ENDIF(HAVE_LIBFFTW)

  IF(QMC_CUDA)
    SET(HAMSRCS_CUDA
      CudaCoulomb.cu 
      NLPP.cu
    )
    SET(HAMSRCS ${HAMSRCS}
      CoulombPBCAA_CUDA.cpp
      CoulombPBCAB_CUDA.cpp
      MPC_CUDA.cpp
      LocalECPotential_CUDA.cpp
      NonLocalECPotential_CUDA.cpp
      CoulombPotential_CUDA.cpp
      SkEstimator_CUDA.cpp
      )
  ENDIF(QMC_CUDA)

ENDIF(OHMMS_DIM MATCHES 3)

IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmcham OBJECT ${HAMSRCS})
  ADD_LIBRARY(qmcham_unit OBJECT ${HAMSRCS})
ELSE()
  ADD_LIBRARY(qmcham ${HAMSRCS})
  ADD_LIBRARY(qmcham_unit ${HAMSRCS})
ENDIF()
IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(qmcham_cuda ${HAMSRCS_CUDA})
  TARGET_LINK_LIBRARIES(qmcham PRIVATE qmcham_cuda)
  TARGET_LINK_LIBRARIES(qmcham_unit PRIVATE qmcham_cuda)
ENDIF(QMC_CUDA)

USE_FAKE_RNG(qmcham_unit)
TARGET_LINK_LIBRARIES(qmcham_unit PUBLIC qmcfakerng)

TARGET_INCLUDE_DIRECTORIES(qmcham PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_INCLUDE_DIRECTORIES(qmcham_unit PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

TARGET_LINK_LIBRARIES(qmcham PUBLIC qmcwfs)
TARGET_LINK_LIBRARIES(qmcham_unit PUBLIC qmcwfs)

TARGET_LINK_LIBRARIES(qmcham PRIVATE einspline Math::BLAS_LAPACK Math::FFTW3)
TARGET_LINK_LIBRARIES(qmcham_unit PRIVATE einspline Math::BLAS_LAPACK Math::FFTW3)
