SET(WFBASE_SRCS
    TrialWaveFunction.cpp
    OrbitalBase.cpp
    OrbitalBuilderBase.cpp
    WaveFunctionFactory.cpp
    ComboOrbital.cpp
    OrbitalConstraintsBase.cpp
   )

SET(JASTROW_SRCS
  Jastrow/LRTwoBodyJastrow.cpp
  Jastrow/JAABuilder.cpp
  Jastrow/JABBuilder.cpp
  Jastrow/PadeJastrowBuilder.cpp
  Jastrow/WMJastrowBuilder.cpp
  Jastrow/AnyConstraints.cpp
  Jastrow/CBSOBuilder.cpp
  Jastrow/JastrowBuilder.cpp
  Jastrow/BsplineJastrowBuilder.cpp
  Jastrow/kSpaceJastrow.cpp
  Jastrow/kSpaceJastrowBuilder.cpp
  Jastrow/RPAJastrow.cpp
  Jastrow/singleRPAJastrowBuilder.cpp
)

#  Jastrow/BsplineFunctor.cpp
# common ferminoic wavefunctions
SET(FERMION_SRCS
  SPOSetBase.cpp
  BasisSetFactory.cpp
  Fermion/DiracDeterminantBase.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  PlaneWave/PWBasis.cpp
  PlaneWave/PWParameterSet.cpp
  PlaneWave/PWOrbitalBuilder.cpp
)

# plane-wave / HEG
IF(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
      ElectronGas/ElectronGasComplexOrbitalBuilder.cpp
      PlaneWave/PWOrbitalSet.cpp
     )
ELSE(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
      ElectronGas/ElectronGasOrbitalBuilder.cpp
      PlaneWave/PWRealOrbitalSet.cpp
     )
ENDIF(QMC_COMPLEX)

# wavefunctions only availbale to 3-dim problems
IF(OHMMS_DIM MATCHES 3)
  SET(FERMION_SRCS ${FERMION_SRCS}
      Jastrow/ThreeBodyGeminal.cpp
      Jastrow/ThreeBodyBlockSparse.cpp
      Jastrow/JastrowBasisBuilder.cpp
      Bspline3DSetBase.cpp
      Bspline3DSet.cpp
      Bspline3DSetTrunc.cpp
      TricubicBsplineSetBuilder.cpp
      TricubicBsplineSetBuilder.1.cpp
      TricubicBsplineSetBuilder.2.cpp
      MolecularOrbitals/STOBuilder.cpp
      MolecularOrbitals/GTOBuilder.cpp
      MolecularOrbitals/NGOBuilder.cpp
      MolecularOrbitals/BsplineAOBuilder.cpp
     )
  IF(HAVE_EINSPLINE)
    SET(FERMION_SRCS 
      ${FERMION_SRCS}
      EinsplineSet.cpp
      EinsplineSetBuilder.cpp
      EinsplineWrapper.cpp
      MuffinTin.cpp
     )
  ENDIF(HAVE_EINSPLINE)
    
  IF(NOT QMC_COMPLEX)
    SET(FERMION_SRCS ${FERMION_SRCS} AGPDeterminant.cpp AGPDeterminantBuilder.cpp)
  ENDIF(NOT QMC_COMPLEX)
ENDIF(OHMMS_DIM MATCHES 3)

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})

#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)
