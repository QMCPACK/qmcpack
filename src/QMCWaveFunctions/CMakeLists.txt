SET(WFSRCS
  OrbitalBuilderBase.cpp
  JAABuilder.cpp
  JABBuilder.cpp
  NJAABuilder.cpp
  NJABBuilder.cpp
  LRTwoBodyJastrow.cpp
  JAAPBCBuilder.cpp
  BasisSetBase.cpp
  ComboOrbital.cpp
  OrbitalConstraintsBase.cpp
  PadeConstraints.cpp
  TwoBodyJastrowBuilder.cpp
  RPAConstraints.cpp
  WMConstraints.cpp
  WaveFunctionFactory.cpp
  SPOSetBase.cpp
  BasisSetFactory.cpp
  Fermion/DiracDeterminantBase.cpp
  Fermion/SlaterDet.cpp
  Fermion/SlaterDetBuilder.cpp
  MolecularOrbitals/STOBuilder.cpp
  MolecularOrbitals/GTOBuilder.cpp
  MolecularOrbitals/NGOBuilder.cpp
)

IF(QMC_COMPLEX)
  SET(WFSRCS ${WFSRCS}
    ElectronGasComplexOrbitalBuilder.cpp
    PlaneWaveOrbitalBuilder.cpp
    PlaneWave/PWOrbitalBuilder.cpp
    PlaneWave/PWOrbitalSet.cpp
  )
ELSE(QMC_COMPLEX)
  SET(WFSRCS ${WFSRCS}
    NumericalOrbitalSetBuilder.cpp
    SplineSetBuilder.cpp
    TricubicBsplineSetBuilder.cpp
    ElectronGasOrbitalBuilder.cpp
    PlaneWaveOrbitalBuilder.cpp
    ThreeBodyGeminal.cpp
    ThreeBodyGeminalBuilder.cpp
    ThreeBodyPadeBuilder.cpp
    AGPDeterminant.cpp
    AGPDeterminantBuilder.cpp
    MolecularOrbitals/STOMolecularOrbitals.cpp
    MolecularOrbitals/GTOMolecularOrbitals.cpp
    MolecularOrbitals/RGFBuilderBase.cpp
    MolecularOrbitals/NumericalRGFBuilder.cpp
    MolecularOrbitals/STO2GridBuilder.cpp
    MolecularOrbitals/GTO2GridBuilder.cpp
    MolecularOrbitals/Any2GridBuilder.cpp
    MolecularOrbitals/GridMolecularOrbitals.cpp
    MolecularOrbitals/MolecularOrbitalBuilder.cpp
    AtomicOrbitals/HFAtomicSTOSet.cpp 
    AtomicOrbitals/HFAtomicSTOSetBuilder.cpp 
    AtomicOrbitals/HeSTOClementiRottie.cpp
  )
ENDIF(QMC_COMPLEX)

####################################
# create libqmcwfs
####################################
ADD_LIBRARY(qmcwfs ${WFSRCS})

