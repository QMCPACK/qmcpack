SET(HAMSRCS
  CoulombPBCAATemp.cpp
  CoulombPBCABTemp.cpp
  LocalCorePolPotential.cpp
  HamiltonianFactory.cpp
  Ylm.cpp
  LocalECPotential.cpp
  NonLocalECPComponent.cpp
  NonLocalECPotential.cpp
  ECPotentialBuilder.cpp
  ECPComponentBuilder.cpp
  NonLocalTOperator.cpp
)

IF(HAVE_LIBFFTW)
  SET(HAMSRCS ${HAMSRCS}
    DispersionRelation.cpp
    MomentumDistribution.cpp
)
ENDIF(HAVE_LIBFFTW)

ADD_LIBRARY(qmcham ${HAMSRCS})
