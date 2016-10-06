#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign  
#//
#// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


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
  QMCFactory/OneDimGridFactory.cpp \
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
  LongRange/KContainer.cpp \
  LongRange/StructFact.cpp \
  LongRange/LPQHIBasis.cpp \
  LongRange/LRCoulombSingleton.cpp\
  LongRange/LRJastrowSingleton.cpp \
  ParticleTags.cpp  \
  ParticleIO/ParticleLayoutIO.cpp  \
  ParticleIO/HDFParticleIO.cpp  \
  ParticleIO/XMLParticleIO.cpp  

 
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
	$(QMCPP_UTILITIES) 

libqmcbase_a_SOURCES = \
	$(QMCPP_PARTICLE)   

INCLUDES = $(OOMPI_INC)

#important: build order,
#  i) current directory for the base libraries
# ii) application-dependent libraries: FortranLib, OOMPI, SlaterKosterTB
#iii) the main application: OhmmsApp
# iv) ready to install the header files and libraries
SUBDIRS = $(oompi) . QMCWaveFunctions QMCHamiltonians QMCDrivers QMCApp
