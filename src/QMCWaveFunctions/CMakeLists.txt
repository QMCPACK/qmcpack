SET(WFBASE_SRCS
  OrbitalBase.cpp
  DiffOrbitalBase.cpp
  OrbitalBuilderBase.cpp
  BasisSetBuilder.cpp
  ProductOrbital.cpp
  OrbitalConstraintsBase.cpp
  SPOInfo.cpp
  SPOSetInfo.cpp
  SPOSetInputInfo.cpp
  SPOSetBase.cpp
  FermionBase.cpp
  OptimizableSPOSet.cpp
  AFMSPOSet.cpp
  CompositeSPOSet.cpp
  HarmonicOscillator/SHOSet.cpp
  HarmonicOscillator/SHOSetBuilder.cpp
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
  Jastrow/JAABuilder.cpp
  Jastrow/JABBuilder.cpp
  IonOrbital.cpp
  IonOrbitalBuilder.cpp
  OptimizableSPOBuilder.cpp
  AFMSPOBuilder.cpp
  Fermion/SPOSetProxy.cpp
  Fermion/SPOSetProxyForMSD.cpp
  )

IF(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasComplexOrbitalBuilder.cpp
    )
ELSE(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasOrbitalBuilder.cpp
    )

ENDIF(QMC_COMPLEX)


# wavefunctions only availbale to 3-dim problems
IF(OHMMS_DIM MATCHES 3)

  SET(JASTROW_SRCS ${JASTROW_SRCS}
    Jastrow/eeI_JastrowBuilder.cpp
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

  IF(QMC_CUDA)
    SET(FERMION_SRCS ${FERMION_SRCS}
      Fermion/determinant_update.cu 
      Fermion/DiracDeterminantCUDA.cpp
      Fermion/SlaterDetCUDA.cpp
      TrialWaveFunction_CUDA.cpp
      )
    SET(JASTROW_SRCS ${JASTROW_SRCS}
      Jastrow/BsplineJastrowCuda.cu
      Jastrow/BsplineJastrowCudaPBC.cu
      Jastrow/OneBodyJastrowOrbitalBspline.cpp
      Jastrow/TwoBodyJastrowOrbitalBspline.cpp
      )
    IF(HAVE_EINSPLINE)
      SET(FERMION_SRCS ${FERMION_SRCS}
        EinsplineSetCuda.cpp
        PhaseFactors.cu
        AtomicOrbitalCuda.cu
        )
    ENDIF(HAVE_EINSPLINE)
  ENDIF(QMC_CUDA)


  IF(HAVE_EINSPLINE)
    SET(FERMION_SRCS ${FERMION_SRCS}
      EinsplineSet.cpp
      EinsplineSetLocal.cpp
      EinsplineSetBuilderCommon.cpp
      EinsplineSetBuilderOld.cpp
      MuffinTin.cpp
      AtomicOrbital.cpp
      EinsplineSetBuilderReadBands_ESHDF.cpp
      EinsplineSetBuilderESHDF.fft.cpp
      EinsplineSetBuilder_createSPOs.cpp
      BandInfo.cpp
      BsplineReaderBase.cpp
      )
    if(NOT IBM_COMPILER)
      SET(FERMION_SRCS ${FERMION_SRCS} EinsplineWrapper.cpp)
    endif(NOT IBM_COMPILER)
  ENDIF(HAVE_EINSPLINE)

  #  IF(QMC_BUILD_LEVEL GREATER 1)
  #
  #    SET(FERMION_SRCS ${FERMION_SRCS}
  #      Bspline3DSetBase.cpp
  #      Bspline3DSet.cpp
  #      Bspline3DSetTrunc.cpp
  #      TricubicBsplineSetBuilder.cpp
  #      TricubicBsplineSetBuilder.1.cpp
  #      TricubicBsplineSetBuilder.2.cpp
  #      PlaneWave/PWBasis.cpp
  #      PlaneWave/PWParameterSet.cpp
  #      PlaneWave/PWOrbitalBuilder.cpp
  #      )
  #  ENDIF(QMC_BUILD_LEVEL GREATER 1)

  #only experimental version
  IF(QMC_BUILD_LEVEL GREATER 2)
   IF(NOT QMC_COMPLEX) 
   SET(FERMION_SRCS ${FERMION_SRCS}
      AGPDeterminant.cpp AGPDeterminantBuilder.cpp
      )
  ENDIF(NOT QMC_COMPLEX)
  ENDIF(QMC_BUILD_LEVEL GREATER 2)

ENDIF(OHMMS_DIM MATCHES 3)

SET(FERMION_SRCS ${FERMION_SRCS}
  Fermion/DiracDeterminantBase.cpp
  Fermion/DiracDeterminantOpt.cpp
  Fermion/DiracDeterminantAFM.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  Fermion/MultiSlaterDeterminant.cpp
  Fermion/MultiSlaterDeterminantFast.cpp
  Fermion/MultiDiracDeterminantBase.cpp
  Fermion/BackflowBuilder.cpp 
  Fermion/DiracDeterminantWithBackflow.cpp
  Fermion/SlaterDetWithBackflow.cpp
  Fermion/MultiSlaterDeterminantWithBackflow.cpp
  BasisSetFactory.cpp
  TrialWaveFunction.cpp
  WaveFunctionFactory.cpp
  )

IF(NOT QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    Fermion/RNDiracDeterminantBase.cpp
    Fermion/RNDiracDeterminantBaseAlternate.cpp
    )
ENDIF(NOT QMC_COMPLEX)

####################################
# create libqmcwfs
####################################
IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
ELSE(QMC_CUDA)
  ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
ENDIF(QMC_CUDA)
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)
