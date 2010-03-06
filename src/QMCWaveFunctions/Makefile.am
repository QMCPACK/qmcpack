#add mpi-related sources and oompi directory for --enable-mpi=yes
WFBASE_SRCS = \
  TrialWaveFunction.cpp \
  OrbitalBuilderBase.cpp \
  WaveFunctionFactory.cpp

JASTROW_SRCS = \
  Jastrow/LRTwoBodyJastrow.cpp \
  Jastrow/JAABuilder.cpp \
  Jastrow/JABBuilder.cpp \
  Jastrow/PadeConstraints.cpp \
  Jastrow/RPAConstraints.cpp \
  Jastrow/AnyConstraints.cpp \
  Jastrow/CBSOBuilder.cpp \
  Jastrow/JastrowBasisBuilder.cpp \
  Jastrow/JastrowBuilder.cpp \
  OrbitalConstraintsBase.cpp \
  ComboOrbital.cpp

FERMION_SRCS =                                  \
  BasisSetBase.cpp                              \
  SPOSetBase.cpp                                \
  BasisSetFactory.cpp                           \
  SplineSetBuilder.cpp                          \
  TricubicBsplineSetBuilder.cpp                 \
  Fermion/DiracDeterminantBase.cpp              \
  Fermion/SlaterDet.cpp                         \
  Fermion/SlaterDetBuilder.cpp                  \
  OptimizableOrbitalSet.cpp                     \
  PlaneWave/PWBasis.cpp                         \
  PlaneWave/PWOrbitalSet.cpp                    \
  PlaneWave/PWParameterSet.cpp                  \
  PlaneWave/PWOrbitalBuilder.cpp                \
  MolecularOrbitals/STOBuilder.cpp              \
  MolecularOrbitals/GTOBuilder.cpp              \
  MolecularOrbitals/NGOBuilder.cpp              \
  MolecularOrbitals/STOMolecularOrbitals.cpp    \
  MolecularOrbitals/GTOMolecularOrbitals.cpp    \
  MolecularOrbitals/RGFBuilderBase.cpp          \
  MolecularOrbitals/NumericalRGFBuilder.cpp     \
  MolecularOrbitals/STO2GridBuilder.cpp         \
  MolecularOrbitals/GTO2GridBuilder.cpp         \
  MolecularOrbitals/Any2GridBuilder.cpp         \
  MolecularOrbitals/GridMolecularOrbitals.cpp   \
  MolecularOrbitals/MolecularOrbitalBuilder.cpp \
  Jastrow/ThreeBodyGeminal.cpp                  \
  Jastrow/ThreeBodyGeminalBuilder.cpp           \
  Jastrow/ThreeBodyPadeBuilder.cpp              \
  AGPDeterminant.cpp                            \
  AGPDeterminantBuilder.cpp

# build the base library used by any application
lib_LIBRARIES = libqmcwfs.a
libqmcwfs_a_SOURCES = \
	$(WFBAS_SRCS) \
	$(JASTROW_SRCS) \
	$(FERMION_SRCS)
