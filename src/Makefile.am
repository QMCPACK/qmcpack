pkginclude_HEADERS= \
	ohmms-config.h 

OHMMS_UTILITIES = \
	Utilities/OhmmsObject.cpp \
	Utilities/OhmmsObject.h \
	Utilities/OhmmsInform.cpp \
	Utilities/OhmmsInform.h \
	Utilities/OhmmsInfo.cpp \
	Utilities/OhmmsInfo.h \
	Utilities/SimpleParser.cpp \
	Utilities/SimpleParser.h \
	Utilities/SpeciesSet.cpp \
	Utilities/SpeciesSet.h \
	Utilities/RandRandom.cpp \
	Utilities/RandRandom.h \
	Utilities/RandomGenerator.cpp \
	Utilities/RandomGenerator.h
	OhmmsApp/ProjectData.cpp \
	OhmmsApp/ProjectData.h \
        OhmmsApp/RandomNumberControl.cpp \
        OhmmsApp/RandomNumberControl.h \
        Platforms/sysutil.cpp\
        Platforms/sysutil.h \
        Optimize/Minimize.cpp\
        Optimize/Minimize.h

   
OHMMS_MESSAGE = \
	Message/Communicate.cpp \
	Message/Communicate.h

OHMMS_PARTICLE = \
	ParticleTags.cpp \
	Particle/ParticleSet.cpp \
	Particle/ParticleSet.h \
	Particle/MCWalkerConfiguration.cpp \
	Particle/MCWalkerConfiguration.h \
        Particle/DistanceTable.cpp \
        Particle/DistanceTable.h \
        Particle/HDFWalkerIO.cpp \
        Particle/HDFWalkerIO.h \
        ParticleIO/ParticleLayoutIO.cpp \
        ParticleIO/XMLParticleIO.h \
        ParticleIO/ParticleLayoutIO.cpp \
        ParticleIO/ParticleLayoutIO.h \
        ParticleIO/HDFParticleIO.cpp \
        ParticleIO/HDFParticleIO.h

# build the base library used by any application
lib_LIBRARIES = libqmcbase.a

libqmcbase_a_SOURCES = \
	$(OHMMS_UTILITIES) \
	$(OHMMS_MESSAGE) \
	$(OHMMS_PARTICLE) 

#important: build order,
#  i) current directory for the base libraries
# ii) application-dependent libraries: FortranLib, OOMPI, SlaterKosterTB
#iii) the main application: OhmmsApp
# iv) ready to install the header files and libraries
SUBDIRS = .  QMC 
