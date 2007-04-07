SET(HAMSRCS
  QMCHamiltonian.cpp
  CoulombPBCAATemp.cpp
  CoulombPBCABTemp.cpp
  LocalCorePolPotential.cpp
  HamiltonianFactory.cpp
  LocalECPotential.cpp
  NonLocalECPComponent.cpp
  NonLocalECPotential.cpp
  ECPotentialBuilder.cpp
  ECPComponentBuilder.cpp
  ECPComponentBuilder.1.cpp
  ECPComponentBuilder.2.cpp
  NonLocalTOperator.cpp
)
#  FSAtomPseudoPot.cpp

IF(HAVE_LIBFFTW)
  SET(HAMSRCS ${HAMSRCS}
    DispersionRelation.cpp
    MomentumDistribution.cpp
)
ENDIF(HAVE_LIBFFTW)

ADD_LIBRARY(qmcham ${HAMSRCS})
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcham STATIC ${HAMSRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcham SHARED ${HAMSRCS})
#ENDIF(QMC_BUILD_STATIC)
