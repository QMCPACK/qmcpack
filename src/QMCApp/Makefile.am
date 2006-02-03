QMC_COMMON_SRCS = \
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp \
  ../QMCWaveFunctions/JAABuilder.cpp \
  ../QMCWaveFunctions/NJAABuilder.cpp \
  ../QMCWaveFunctions/JABBuilder.cpp \
  ../QMCWaveFunctions/NJABBuilder.cpp \
  ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.cpp \
  ../QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.cpp \
  ../QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp\
  ../QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp\
  ../QMCWaveFunctions/ElectronGasOrbitalBuilder.cpp\
  ../QMCWaveFunctions/ThreeBodyGeminal.cpp\
  ../QMCWaveFunctions/ThreeBodyGeminalBuilder.cpp\
  ../QMCHamiltonians/LocalPPotential.cpp\
  ../QMCHamiltonians/NonLocalPPotential.cpp\
  ../QMCHamiltonians/LocalCorePolPotential.cpp \
  ../QMCWaveFunctions/WaveFunctionFactory.cpp \
  ../QMCHamiltonians/HamiltonianFactory.cpp

QMC_DRIVERS = \
  ../QMCDrivers/SimpleFixedNodeBranch.cpp\
  ../QMCDrivers/QMCDriver.cpp\
  ../QMCDrivers/DummyQMC.cpp\
  ../QMCDrivers/VMC.cpp\
  ../QMCDrivers/VMCParticleByParticle.cpp\
  ../QMCDrivers/MolecuDMC.cpp\
  ../QMCDrivers/DMCParticleByParticle.cpp\
  ../QMCDrivers/WalkerControlBase.cpp\
  ../QMCDrivers/QMCOptimize.cpp\
  ../QMCDrivers/QMCCostFunction.cpp\
  ../QMCDrivers/DMCPbyPOpenMP.cpp \
  ../QMCDrivers/DMCPbyPUpdate.cpp \
  ../QMCDrivers/DMCPbyP.cpp \
  ../Estimators/ScalarEstimatorManager.cpp

QMC_APP =\
  ParticleSetPool.cpp\
  WaveFunctionPool.cpp\
  HamiltonianPool.cpp\
  InitMolecularSystem.cpp\
  QMCAppBase.cpp\
  QMCMain.cpp\
  qmcapp.cpp
if NEEDTOBUILDOOMPI
QMC_MPP_SRC = \
  ../QMCDrivers/GlobalWalkerControl.cpp
endif

bin_PROGRAMS = qmcapp
qmcapp_SOURCES = \
$(QMC_COMMON_SRCS) \
$(QMC_DRIVERS) \
$(QMC_MPP_SRC) \
$(QMC_APP)

if NEEDTOBUILDOOMPI
OOMPI_INC = -I$(top_srcdir)/src/OOMPI
OOMPI_LIB = ${top_builddir}/src/OOMPI/liboompi.a
endif

INCLUDES = -I$(top_srcdir)/src $(OOMPI_INC)

qmcapp_LDADD = @QMCBASE_LIBS@  $(OOMPI_LIB)
qmcapp_LDFLAGS = @LDFLAGS@ 
qmcapp_DEPENDENCIES = ${top_builddir}/src/libqmcbase.a $(OOMPI_LIB)
