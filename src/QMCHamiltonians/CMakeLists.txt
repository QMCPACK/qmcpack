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
  QMCHamiltonianBase.cpp
  QMCHamiltonian.cpp
  CoulombPBCAA.cpp
  CoulombPBCAB.cpp
  NonLocalTOperator.cpp
  ForwardWalking.cpp
  PairCorrEstimator.cpp
  LocalMomentEstimator.cpp
  DensityEstimator.cpp
  SkPot.cpp
  SkEstimator.cpp
  SkAllEstimator.cpp
  MomentumEstimator.cpp
  ForceBase.cpp
  HamiltonianFactory.cpp
  CoulombPotentialFactory.cpp
  HarmonicExternalPotential.cpp
  StaticStructureFactor.cpp
  OrbitalImages.cpp
  SpinDensity.cpp
  SpeciesKineticEnergy.cpp
  LatticeDeviationEstimator.cpp
  )

IF(OHMMS_DIM MATCHES 3)
  SET(HAMSRCS ${HAMSRCS}
    ChiesaCorrection.cpp
    PulayForce.cpp
    ZeroVarianceForce.cpp
    ForceCeperley.cpp
    ForceChiesaPBCAA.cpp
#    StressPBCAA.cpp
#    StressPBCAB.cpp
#    StressKinetic.cpp
    StressPBC.cpp
    LocalCorePolPotential.cpp
    LocalECPotential.cpp
    NonLocalECPComponent.cpp
    NonLocalECPotential.deriv.cpp
    NonLocalECPotential.cpp
    ECPotentialBuilder.cpp
    ECPComponentBuilder.cpp
    ECPComponentBuilder.1.cpp
    ECPComponentBuilder.2.cpp
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
    SET(HAMSRCS ${HAMSRCS}
      CudaCoulomb.cu 
      NLPP.cu
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

IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(qmcham ${HAMSRCS})
ELSE(QMC_CUDA)
  ADD_LIBRARY(qmcham ${HAMSRCS})
ENDIF(QMC_CUDA)


#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcham STATIC ${HAMSRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcham SHARED ${HAMSRCS})
#ENDIF(QMC_BUILD_STATIC)
