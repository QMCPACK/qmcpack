SET(HAMSRCS
  QMCHamiltonianBase.cpp
  QMCHamiltonian.cpp
  CoulombPBCAA.cpp
  CoulombPBCAB.cpp
  NonLocalTOperator.cpp
  NumericalRadialPotential.cpp
  ForwardWalking.cpp
  PairCorrEstimator.cpp
  LocalMomentEstimator.cpp
  DensityEstimator.cpp
  SkPot.cpp
  SkEstimator.cpp
  MomentumEstimator.cpp
  ForceBase.cpp
  HamiltonianFactory.cpp
  CoulombPotentialFactory.cpp
  ../Estimators/TraceManager.cpp
  HarmonicExternalPotential.cpp
  SpinDensity.cpp
  StaticStructureFactor.cpp
  ../Estimators/PostProcessor.cpp
  ../Estimators/SpinDensityPostProcessor.cpp
  OrbitalImages.cpp
  DensityMatrices1B.cpp
)
#  FSAtomPseudoPot.cpp

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
    NonLocalECPComponent.vp.cpp
    NonLocalECPotential.deriv.cpp
    NonLocalECPotential.cpp
    ECPotentialBuilder.cpp
    ECPComponentBuilder.cpp
    ECPComponentBuilder.1.cpp
    ECPComponentBuilder.2.cpp
    ReferencePoints.cpp
    SpaceGrid.cpp
    EnergyDensityEstimator.cpp
    )

  IF(HAVE_LIBFFTW)
    SET(HAMSRCS ${HAMSRCS} MPC.cpp VHXC.cpp)
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

  # Add research-related sources here
  IF(QMC_BUILD_LEVEL GREATER 1)
    SET(HAMSRCS ${HAMSRCS}
      LennardJones_smoothed.cpp
      HFDHE2_Moroni1995.cpp
      HeSAPT_smoothed.cpp
     )
  ENDIF(QMC_BUILD_LEVEL GREATER 1)
ENDIF(OHMMS_DIM MATCHES 3)

#IF(HAVE_LIBFFTW)
#  SET(HAMSRCS ${HAMSRCS}
#    DispersionRelation.cpp
#    MomentumDistribution.cpp
#)
#ENDIF(HAVE_LIBFFTW)

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
