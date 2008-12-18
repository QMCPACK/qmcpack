SET(HAMSRCS
  QMCHamiltonian.cpp
  CoulombPBCAATemp.cpp
  CoulombPBCABTemp.cpp
  HamiltonianFactory.cpp
  NonLocalTOperator.cpp
  NumericalRadialPotential.cpp
  RPAPressure.cpp
  ForceBase.cpp
  ForceCeperley.cpp
  ChiesaCorrection.cpp
  PairCorrEstimator.cpp
  DensityEstimator.cpp
)
#  FSAtomPseudoPot.cpp

IF(OHMMS_DIM MATCHES 3)
  SET(HAMSRCS ${HAMSRCS}
      LocalCorePolPotential.cpp
      LocalECPotential.cpp
      NonLocalECPComponent.cpp
      NonLocalECPotential.cpp
      ECPotentialBuilder.cpp
      ECPComponentBuilder.cpp
      ECPComponentBuilder.1.cpp
      ECPComponentBuilder.2.cpp
     )
ENDIF(OHMMS_DIM MATCHES 3)

#IF(HAVE_LIBFFTW)
#  SET(HAMSRCS ${HAMSRCS}
#    DispersionRelation.cpp
#    MomentumDistribution.cpp
#)
#ENDIF(HAVE_LIBFFTW)

IF(HAVE_LIBFFTW)
  SET(HAMSRCS ${HAMSRCS} MPC.cpp)
ENDIF(HAVE_LIBFFTW)

ADD_LIBRARY(qmcham ${HAMSRCS})
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcham STATIC ${HAMSRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcham SHARED ${HAMSRCS})
#ENDIF(QMC_BUILD_STATIC)
