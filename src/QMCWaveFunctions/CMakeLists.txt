
SET(WFBASE_SRCS
  TrialWaveFunction.cpp
  OrbitalBuilderBase.cpp
  WaveFunctionFactory.cpp
  )

SET(JASTROW_SRCS
  JAABuilder.cpp
  JABBuilder.cpp
  NJAABuilder.cpp
  NJABBuilder.cpp
  LRTwoBodyJastrow.cpp
  JAAPBCBuilder.cpp
  TwoBodyJastrowBuilder.cpp
  OrbitalConstraintsBase.cpp
  PadeConstraints.cpp
  RPAConstraints.cpp
  WMConstraints.cpp
  ComboOrbital.cpp
)

SET(FERMION_SRCS
  BasisSetBase.cpp
  SPOSetBase.cpp
  BasisSetFactory.cpp
  SplineSetBuilder.cpp
  TricubicBsplineSetBuilder.cpp
  Fermion/DiracDeterminantBase.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  PlaneWave/PWOrbitalBuilder.cpp
  PlaneWave/PWOrbitalSet.cpp
  MolecularOrbitals/STOBuilder.cpp
  MolecularOrbitals/GTOBuilder.cpp
  MolecularOrbitals/NGOBuilder.cpp
  MolecularOrbitals/STOMolecularOrbitals.cpp
  MolecularOrbitals/GTOMolecularOrbitals.cpp
  MolecularOrbitals/RGFBuilderBase.cpp
  MolecularOrbitals/NumericalRGFBuilder.cpp
  MolecularOrbitals/STO2GridBuilder.cpp
  MolecularOrbitals/GTO2GridBuilder.cpp
  MolecularOrbitals/Any2GridBuilder.cpp
  MolecularOrbitals/GridMolecularOrbitals.cpp
  MolecularOrbitals/MolecularOrbitalBuilder.cpp
)

IF(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGasComplexOrbitalBuilder.cpp
  )
ELSE(QMC_COMPLEX)
  SET(FERMION_SRCS ${FERMION_SRCS}
    ElectronGasOrbitalBuilder.cpp
    ThreeBodyGeminal.cpp
    ThreeBodyGeminalBuilder.cpp
    ThreeBodyPadeBuilder.cpp
    AGPDeterminant.cpp
    AGPDeterminantBuilder.cpp
  )
ENDIF(QMC_COMPLEX)

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${WFBASE_SRCS} ${JASTROW_SRCS} ${FERMION_SRCS})

