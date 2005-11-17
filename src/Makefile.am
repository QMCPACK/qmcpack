pkginclude_HEADERS= \
	ohmms-config.h 

QMCPP_UTILITIES = \
	Utilities/OhmmsObject.cpp \
	Utilities/OhmmsInform.cpp \
	Utilities/OhmmsInfo.cpp \
	Utilities/SimpleParser.cpp \
	Utilities/SpeciesSet.cpp \
	Utilities/RandRandom.cpp \
	Utilities/RandomGenerator.cpp \
	OhmmsData/Libxml2Doc.cpp \
	OhmmsApp/ProjectData.cpp \
        OhmmsApp/RandomNumberControl.cpp \
        Platforms/sysutil.cpp\
	Numerics/OhmmsBlas.cpp

   
QMCPP_MESSAGE = \
	Message/Communicate.cpp 

QMCPP_PARTICLE = \
	ParticleTags.cpp \
        Lattice/Uniform3DGridLayout.cpp \
        Particle/ParticleSet.BC.cpp \
	Particle/ParticleSet.cpp \
	Particle/MCWalkerConfiguration.cpp \
        Particle/DistanceTable.cpp \
        Particle/HDFWalkerInput0.cpp \
        Particle/HDFWalkerInputCollect.cpp \
        Particle/HDFWalkerInputManager.cpp \
        Particle/HDFWalkerOutput.cpp \
        ParticleIO/XMLParticleIO.cpp \
        ParticleIO/ParticleLayoutIO.cpp \
        ParticleIO/HDFParticleIO.cpp \
        LongRange/KContainer.cpp \
        LongRange/StructFact.cpp 

QMC_BASE = \
	QMCHamiltonians/QMCHamiltonian.cpp\
	QMCWaveFunctions/JastrowBuilder.cpp\
	QMCWaveFunctions/TrialWaveFunction.cpp

#add mpi-related sources and oompi directory for --enable-mpi=yes
if NEEDTOBUILDOOMPI
oompi=OOMPI
OOMPI_INC =  -I$(top_srcdir)/src/OOMPI
endif

INCLUDES = $(OOMPI_INC)

# build the base library used by any application
lib_LIBRARIES = libqmcbase.a

libqmcbase_a_SOURCES = \
	$(QMCPP_UTILITIES) \
	$(QMCPP_MESSAGE) \
	$(QMCPP_PARTICLE) \
	$(QMC_BASE)

#important: build order,
#  i) current directory for the base libraries
# ii) application-dependent libraries: FortranLib, OOMPI, SlaterKosterTB
#iii) the main application: OhmmsApp
# iv) ready to install the header files and libraries
SUBDIRS = .  $(oompi) QMCApp
