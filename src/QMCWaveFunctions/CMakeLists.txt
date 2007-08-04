SET(WFBASE_SRCS
    TrialWaveFunction.cpp
    OrbitalBuilderBase.cpp
    WaveFunctionFactory.cpp
    ComboOrbital.cpp
    OrbitalConstraintsBase.cpp
   )

SET(JASTROW_SRCS
  Jastrow/LRTwoBodyJastrow.cpp
  Jastrow/JAABuilder.cpp
  Jastrow/JABBuilder.cpp
  Jastrow/PadeConstraints.cpp
  Jastrow/RPAConstraints.cpp
  Jastrow/AnyConstraints.cpp
  Jastrow/PolyConstraints.cpp
  Jastrow/CBSOBuilder.cpp
  Jastrow/JastrowBasisBuilder.cpp
  Jastrow/JastrowBuilder.cpp
  Jastrow/ThreeBodyGeminal.cpp
  Jastrow/ThreeBodyBlockSparse.cpp
)

SET(FERMION_SRCS
  SPOSetBase.cpp
  BasisSetFactory.cpp
  Bspline3DSetBase.cpp
  Bspline3DSet.cpp
  Bspline3DSetTrunc.cpp
  TricubicBsplineSetBuilder.cpp
  TricubicBsplineSetBuilder.1.cpp
  TricubicBsplineSetBuilder.2.cpp
  Fermion/DiracDeterminantBase.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  PlaneWave/PWBasis.cpp
  PlaneWave/PWParameterSet.cpp
  PlaneWave/PWOrbitalBuilder.cpp
  MolecularOrbitals/STOBuilder.cpp
  MolecularOrbitals/GTOBuilder.cpp
  MolecularOrbitals/NGOBuilder.cpp
)

IF(HAVE_EINSPLINE)
  SET(FERMION_SRCS 
      ${FERMION_SRCS}
      EinsplineSet.cpp
  )
ENDIF(HAVE_EINSPLINE)
    
IF(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasComplexOrbitalBuilder.cpp
    PlaneWave/PWOrbitalSet.cpp
  )
ELSE(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasOrbitalBuilder.cpp
    AGPDeterminant.cpp
    AGPDeterminantBuilder.cpp
    PlaneWave/PWRealOrbitalSet.cpp
  )
ENDIF(QMC_COMPLEX)

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs STATIC ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcwfs SHARED ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})
#ENDIF(QMC_BUILD_STATIC)

## THESE ARE NOT USED ANYMORE
#  Jastrow/ThreeBodyGeminalBuilder.cpp
#  Jastrow/ThreeBodyPadeBuilder.cpp
#  MolecularOrbitals/GridMolecularOrbitals.cpp
#  MolecularOrbitals/RGFBuilderBase.cpp
#  MolecularOrbitals/NumericalRGFBuilder.cpp
#  MolecularOrbitals/Any2GridBuilder.cpp
#  MolecularOrbitals/STOMolecularOrbitals.cpp
#  MolecularOrbitals/GTOMolecularOrbitals.cpp
#  MolecularOrbitals/STO2GridBuilder.cpp
#  MolecularOrbitals/GTO2GridBuilder.cpp
#  MolecularOrbitals/MolecularOrbitalBuilder.cpp

