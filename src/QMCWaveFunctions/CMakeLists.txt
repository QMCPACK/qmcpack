SET(WFBASE_SRCS
    TrialWaveFunction.cpp
    OrbitalBase.cpp
    DiffOrbitalBase.cpp
    OrbitalBuilderBase.cpp
    WaveFunctionFactory.cpp
    ProductOrbital.cpp
    OrbitalConstraintsBase.cpp
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
  Jastrow/kSpaceJastrow.cpp
  Jastrow/kSpaceJastrowBuilder.cpp
  Jastrow/RPAJastrow.cpp
  Jastrow/singleRPAJastrowBuilder.cpp
)

#common jastrows for three D
IF(OHMMS_DIM MATCHES 3)
  SET(JASTROW_SRCS ${JASTROW_SRCS}
    Jastrow/ThreeBodyGeminal.cpp
    Jastrow/ThreeBodyBlockSparse.cpp
    Jastrow/JastrowBasisBuilder.cpp
    Jastrow/CBSOBuilder.cpp
    )
ENDIF(OHMMS_DIM MATCHES 3)

#specialized jastrows: extremely taxing for some compilers
IF(QMC_BUILD_COMPLETE)
  SET(JASTROW_SRCS ${JASTROW_SRCS}
    Jastrow/JAABuilder.cpp
    Jastrow/JABBuilder.cpp
    Jastrow/WMJastrowBuilder.cpp
    )
ENDIF(QMC_BUILD_COMPLETE)

# common ferminoic wavefunctions
SET(FERMION_SRCS
  SPOSetBase.cpp
  BasisSetFactory.cpp
  Fermion/DiracDeterminantBase.cpp
  Fermion/DiracDeterminantIterative.cpp
  Fermion/DiracDeterminantTruncation.cpp
  Fermion/MultiDiracDeterminantBase.cpp
  Fermion/ILUGMRESInterface.c
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
)

# wavefunctions only availbale to 3-dim problems
IF(OHMMS_DIM MATCHES 3)
  SET(FERMION_SRCS ${FERMION_SRCS}
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
      AtomicOrbital.cpp
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
    
ENDIF(OHMMS_DIM MATCHES 3)

IF(QMC_BUILD_COMPLETE)
  SET(FERMION_SRCS ${FERMION_SRCS}
    PlaneWave/PWBasis.cpp
    PlaneWave/PWParameterSet.cpp
    PlaneWave/PWOrbitalBuilder.cpp
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
      AGPDeterminant.cpp AGPDeterminantBuilder.cpp
      )
  ENDIF(QMC_COMPLEX)
ENDIF(QMC_BUILD_COMPLETE)

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})

#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)
