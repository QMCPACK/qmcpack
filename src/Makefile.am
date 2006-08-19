#add mpi-related sources and oompi directory for --enable-mpi=yes
pkginclude_HEADERS= \
	ohmms-config.h 

QMCPP_UTILITIES = \
  Utilities/RandRandom.cpp  \
  Utilities/RandomGenerator.cpp  \
  Utilities/OhmmsObject.cpp  \
  Utilities/SpeciesSet.cpp  \
  Utilities/SimpleParser.cpp  \
  Utilities/OhmmsInform.cpp  \
  Utilities/OhmmsInfo.cpp  \
  OhmmsData/Libxml2Doc.cpp \
  OhmmsApp/ProjectData.cpp \
  OhmmsApp/RandomNumberControl.cpp\
  Numerics/OhmmsBlas.cpp \
	Numerics/FFTEngines.h \
  Platforms/sysutil.cpp \
  QMCFactory/OneDimGridFactory.cpp
   
QMCPP_MESSAGE = \
  Message/Communicate.cpp 

QMCPP_PARTICLE = \
  Lattice/Uniform3DGridLayout.cpp  \
  Particle/ParticleSet.cpp  \
  Particle/ParticleSet.BC.cpp  \
  Particle/MCWalkerConfiguration.cpp  \
  Particle/DistanceTable.cpp \
  Particle/HDFWalkerOutput.cpp \
  Particle/HDFWalkerInput0.cpp \
  Particle/HDFWalkerInputCollect.cpp \
  Particle/HDFWalkerInputManager.cpp \
  LongRange/KContainer.cpp \
  LongRange/StructFact.cpp \
  ParticleTags.cpp  \
  ParticleIO/ParticleLayoutIO.cpp  \
  ParticleIO/HDFParticleIO.cpp  \
  ParticleIO/XMLParticleIO.cpp  

QMC_WFS = \
  QMCWaveFunctions/OrbitalBuilderBase.cpp \
  QMCWaveFunctions/JastrowBuilder.cpp\
  QMCWaveFunctions/TrialWaveFunction.cpp  \
  QMCWaveFunctions/JAABuilder.cpp \
  QMCWaveFunctions/NJAABuilder.cpp \
  QMCWaveFunctions/JABBuilder.cpp \
  QMCWaveFunctions/NJABBuilder.cpp \
  QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSet.cpp \
  QMCWaveFunctions/AtomicOrbitals/HFAtomicSTOSetBuilder.cpp \
  QMCWaveFunctions/AtomicOrbitals/HeSTOClementiRottie.cpp\
  QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.cpp\
  QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.cpp\
  QMCWaveFunctions/MolecularOrbitals/RGFBuilderBase.cpp\
  QMCWaveFunctions/MolecularOrbitals/NumericalRGFBuilder.cpp\
  QMCWaveFunctions/MolecularOrbitals/STO2GridBuilder.cpp\
  QMCWaveFunctions/MolecularOrbitals/GTO2GridBuilder.cpp\
  QMCWaveFunctions/MolecularOrbitals/Any2GridBuilder.cpp\
  QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.cpp\
  QMCWaveFunctions/MolecularOrbitals/MolecularOrbitalBuilder.cpp\
  QMCWaveFunctions/ElectronGasOrbitalBuilder.cpp\
  QMCWaveFunctions/ThreeBodyGeminal.cpp\
  QMCWaveFunctions/ThreeBodyGeminalBuilder.cpp\
  QMCWaveFunctions/AGPDeterminant.cpp \
  QMCWaveFunctions/AGPDeterminantBuilder.cpp \
  QMCWaveFunctions/PlaneWaveOrbitalBuilder.cpp \
  QMCWaveFunctions/WaveFunctionFactory.cpp 
 
QMC_HAM = \
  QMCHamiltonians/QMCHamiltonian.cpp\
  QMCHamiltonians/LocalPPotential.cpp\
  QMCHamiltonians/NonLocalPPotential.cpp\
  QMCHamiltonians/LocalCorePolPotential.cpp \
  QMCHamiltonians/LocalECPotential.cpp \
  QMCHamiltonians/NonLocalECPotential.cpp \
  QMCHamiltonians/ECPComponentBuilder.cpp \
  QMCHamiltonians/ECPotentialBuilder.cpp \
  QMCHamiltonians/HamiltonianFactory.cpp

QMC_DRIVERS = \
  QMCDrivers/SimpleFixedNodeBranch.cpp\
  QMCDrivers/QMCDriver.cpp\
  QMCDrivers/DummyQMC.cpp\
  QMCDrivers/VMC.cpp\
  QMCDrivers/VMCParticleByParticle.cpp\
  QMCDrivers/QMCOptimize.cpp\
  QMCDrivers/QMCCostFunction.cpp\
  QMCDrivers/WaveFunctionTester.cpp \
  Estimators/ScalarEstimatorManager.cpp \
  QMCDrivers/WalkerControlBase.cpp \
  QMCDrivers/DMC/DMCMoveAll.cpp \
  QMCDrivers/DMC/DMCPbyP.cpp \
  QMCDrivers/DMC/DMCFactory.cpp \
  QMCDrivers/DMC/DMCUpdateBase.cpp \
  QMCDrivers/DMC/DMCUpdateAll.cpp \
  QMCDrivers/DMC/DMCUpdatePbyP.cpp \
  QMCDrivers/DMC/WalkerControlFactory.cpp \
  QMCDrivers/DMC/WalkerReconfiguration.cpp \
  QMCDrivers/VMCMultiple.cpp \
  QMCDrivers/VMCPbyPMultiple.cpp \
  QMCDrivers/PolymerEstimator.cpp \
  QMCDrivers/RQMCMultiple.cpp \
  QMCDrivers/MultiChain.cpp \
  Estimators/MultipleEnergyEstimator.cpp \
  QMCDrivers/SpaceWarp.cpp \
  QMCDrivers/VMCMultipleWarp.cpp \
  QMCDrivers/VMCPbyPMultiWarp.cpp \
  QMCDrivers/RQMCMultiWarp.cpp

if NEEDTOBUILDOOMPI
oompi=OOMPI
OOMPI_INC =  -I$(top_srcdir)/src/OOMPI
QMC_MPP_SRC = \
  QMCDrivers/DMC/WalkerControlMPI.cpp \
  QMCDrivers/DMC/WalkerReconfigurationMPI.cpp
endif

# build the base library used by any application
lib_LIBRARIES = libqmcutil.a libqmcbase.a 
libqmcutil_a_SOURCES = \
	$(QMCPP_UTILITIES) \
	$(QMCPP_MESSAGE) 

libqmcbase_a_SOURCES = \
	$(QMCPP_PARTICLE)   \
	$(QMC_WFS) \
	$(QMC_HAM) \
	$(QMC_DRIVERS) 

INCLUDES = $(OOMPI_INC)

#important: build order,
#  i) current directory for the base libraries
# ii) application-dependent libraries: FortranLib, OOMPI, SlaterKosterTB
#iii) the main application: OhmmsApp
# iv) ready to install the header files and libraries
SUBDIRS = $(oompi) . QMCApp
