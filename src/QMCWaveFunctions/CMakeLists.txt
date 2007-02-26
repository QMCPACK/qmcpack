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
  Jastrow/CBSOBuilder.cpp
  Jastrow/JastrowBasisBuilder.cpp
  Jastrow/JastrowBuilder.cpp
  Jastrow/ThreeBodyGeminal.cpp
  Jastrow/ThreeBodyBlockSparse.cpp
)

SET(FERMION_SRCS
  SPOSetBase.cpp
  BasisSetFactory.cpp
  SplineSetBuilder.cpp
  TricubicBsplineSetBuilder.cpp
  Fermion/DiracDeterminantBase.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  PlaneWave/PWBasis.cpp
  PlaneWave/PWOrbitalSet.cpp
  PlaneWave/PWParameterSet.cpp
  PlaneWave/PWOrbitalBuilder.cpp
  MolecularOrbitals/STOBuilder.cpp
  MolecularOrbitals/GTOBuilder.cpp
  MolecularOrbitals/NGOBuilder.cpp
)

IF(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasComplexOrbitalBuilder.cpp
  )
ELSE(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGas/ElectronGasOrbitalBuilder.cpp
    AGPDeterminant.cpp
    AGPDeterminantBuilder.cpp
  )
ENDIF(QMC_COMPLEX)

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})

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

