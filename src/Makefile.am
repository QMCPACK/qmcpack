pkginclude_HEADERS= \
	ohmms-config.h 

OHMMS_UTILITIES = \
	Utilities/OhmmsObject.cpp \
	Utilities/OhmmsInform.cpp \
	Utilities/OhmmsInfo.cpp \
	Utilities/SimpleParser.cpp \
	Utilities/SpeciesSet.cpp \
	Utilities/RandRandom.cpp \
	Utilities/RandomGenerator.cpp \
	OhmmsApp/ProjectData.cpp \
        OhmmsApp/RandomNumberControl.cpp \
        Platforms/sysutil.cpp\
        Optimize/Minimize.cpp\
	Numerics/OhmmsBlas.cpp

   
OHMMS_MESSAGE = \
	Message/Communicate.cpp 

OHMMS_PARTICLE = \
	ParticleTags.cpp \
        Lattice/Uniform3DGridLayout.cpp \
        LongRange/KContainer.cpp \
        LongRange/StructFact.cpp \
        Particle/ParticleSet.BC.cpp \
	Particle/ParticleSet.cpp \
	Particle/MCWalkerConfiguration.cpp \
        Particle/DistanceTable.cpp \
        Particle/HDFWalkerIO.cpp \
        ParticleIO/XMLParticleIO.cpp \
        ParticleIO/ParticleLayoutIO.cpp \
        ParticleIO/HDFParticleIO.cpp \
        LongRange/KContainer.cpp \
        LongRange/StructFact.cpp 

QMC_BASE = \
	QMCHamiltonians/QMCHamiltonian.cpp\
	QMCWaveFunctions/JastrowBuilder.cpp\
	QMCWaveFunctions/TrialWaveFunction.cpp

# build the base library used by any application
lib_LIBRARIES = libqmcbase.a

libqmcbase_a_SOURCES = \
	$(OHMMS_UTILITIES) \
	$(OHMMS_MESSAGE) \
	$(OHMMS_PARTICLE) \
	$(QMC_BASE)

#important: build order,
#  i) current directory for the base libraries
# ii) application-dependent libraries: FortranLib, OOMPI, SlaterKosterTB
#iii) the main application: OhmmsApp
# iv) ready to install the header files and libraries
SUBDIRS = .  QMCApp
