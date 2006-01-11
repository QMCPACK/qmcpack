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

QMC_BASE = \
	QMCHamiltonians/QMCHamiltonian.cpp\
	QMCWaveFunctions/JastrowBuilder.cpp\
	QMCWaveFunctions/TrialWaveFunction.cpp

if NEEDTOBUILDOOMPI
oompi=OOMPI
OOMPI_INC =  -I$(top_srcdir)/src/OOMPI
endif

# build the base library used by any application
lib_LIBRARIES = libqmcbase.a


libqmcbase_a_SOURCES = \
	$(QMCPP_UTILITIES) \
	$(QMCPP_MESSAGE) \
	$(QMCPP_PARTICLE) \
	$(QMC_BASE)

INCLUDES = $(OOMPI_INC)

#important: build order,
#  i) current directory for the base libraries
# ii) application-dependent libraries: FortranLib, OOMPI, SlaterKosterTB
#iii) the main application: OhmmsApp
# iv) ready to install the header files and libraries
SUBDIRS = $(oompi) . QMCApp
