SET(WFBASE_SRCS
    OrbitalBase.cpp
    DiffOrbitalBase.cpp
    OrbitalBuilderBase.cpp
    ProductOrbital.cpp
    OrbitalConstraintsBase.cpp
    SPOSetBase.cpp
   )

########################
# build jastrows
########################
#common jastrows
SET(JASTROW_SRCS
  Jastrow/LRTwoBodyJastrow.cpp
  Jastrow/PadeJastrowBuilder.cpp
  Jastrow/JastrowBuilder.cpp
  Jastrow/BsplineJastrowBuilder.cpp
  Jastrow/eeI_JastrowBuilder.cpp
  Jastrow/kSpaceJastrow.cpp
  Jastrow/kSpaceJastrowBuilder.cpp
  Jastrow/RPAJastrow.cpp
  Jastrow/singleRPAJastrowBuilder.cpp
  IonOrbital.cpp
  IonOrbitalBuilder.cpp
)

#specialized jastrows: extremely taxing for some compilers
IF(QMC_BUILD_LEVEL GREATER 1)
  SET(JASTROW_SRCS ${JASTROW_SRCS}
    Jastrow/JAABuilder.cpp
    Jastrow/JABBuilder.cpp
    Jastrow/WMJastrowBuilder.cpp
    )
ENDIF(QMC_BUILD_LEVEL GREATER 1)

# wavefunctions only availbale to 3-dim problems
IF(OHMMS_DIM MATCHES 3)

  SET(JASTROW_SRCS ${JASTROW_SRCS}
    Jastrow/ThreeBodyGeminal.cpp
    Jastrow/ThreeBodyBlockSparse.cpp
    Jastrow/JastrowBasisBuilder.cpp
    Jastrow/CBSOBuilder.cpp
    )

  SET(FERMION_SRCS ${FERMION_SRCS}
    MolecularOrbitals/STOBuilder.cpp
    MolecularOrbitals/GTOBuilder.cpp
    MolecularOrbitals/NGOBuilder.cpp
    MolecularOrbitals/BsplineAOBuilder.cpp
    )

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

  IF(HAVE_EINSPLINE)
    SET(FERMION_SRCS ${FERMION_SRCS}
      EinsplineSet.cpp
      EinsplineSetBuilder.cpp
      EinsplineWrapper.cpp
      MuffinTin.cpp
      AtomicOrbital.cpp
      )
  ENDIF(HAVE_EINSPLINE)


  IF(QMC_BUILD_LEVEL GREATER 1)

    SET(FERMION_SRCS ${FERMION_SRCS}
      Bspline3DSetBase.cpp
      Bspline3DSet.cpp
      Bspline3DSetTrunc.cpp
      TricubicBsplineSetBuilder.cpp
      TricubicBsplineSetBuilder.1.cpp
      TricubicBsplineSetBuilder.2.cpp
      PlaneWave/PWBasis.cpp
      PlaneWave/PWParameterSet.cpp
      PlaneWave/PWOrbitalBuilder.cpp
      )
  ENDIF(QMC_BUILD_LEVEL GREATER 1)

  #only experimental version
  IF(QMC_BUILD_LEVEL GREATER 2)
    SET(FERMION_SRCS ${FERMION_SRCS}
      AGPDeterminant.cpp AGPDeterminantBuilder.cpp
      )
  ENDIF(QMC_BUILD_LEVEL GREATER 2)
    
ENDIF(OHMMS_DIM MATCHES 3)

SET(FERMION_SRCS ${FERMION_SRCS}
  Fermion/DiracDeterminantBase.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  MultiSlaterDeterminant.cpp
  BasisSetFactory.cpp
  TrialWaveFunction.cpp
  WaveFunctionFactory.cpp
  )

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${FERMION_SRCS} ${JASTROW_SRCS} ${WFBASE_SRCS})

#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)
