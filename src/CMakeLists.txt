#CONFIGURE_FILE(${qmcPlusPlus_SOURCE_DIR}/src/Configuration.h.in
#               ${qmcPlusPlus_BINARY_DIR}/src/Configuration.h)
CONFIGURE_FILE(${qmcPlusPlus_SOURCE_DIR}/src/ohmms-config.h.cmake.in
               ${qmcPlusPlus_BINARY_DIR}/src/ohmms-config.h)

SET(UTILITIES
  Utilities/RandRandom.cpp 
  Utilities/RandomGenerator.cpp 
  Utilities/OhmmsObject.cpp 
  Utilities/SpeciesSet.cpp 
  Utilities/SimpleParser.cpp 
  Utilities/OhmmsInform.cpp 
  Utilities/OhmmsInfo.cpp 
  OhmmsApp/ProjectData.cpp
  OhmmsApp/RandomNumberControl.cpp
  Numerics/OhmmsBlas.cpp
  Platforms/sysutil.cpp
)

SET(PARTICLE
  Lattice/Uniform3DGridLayout.cpp 
  Particle/ParticleSet.cpp 
  Particle/ParticleSet.BC.cpp 
  Particle/MCWalkerConfiguration.cpp 
  Particle/DistanceTable.cpp
  Particle/HDFWalkerIO.cpp
  LongRange/KContainer.cpp
  LongRange/StructFact.cpp
)

SET(PARTICLEIO
  ParticleTags.cpp 
  ParticleIO/ParticleLayoutIO.cpp 
  ParticleIO/XMLParticleIO.cpp 
  ParticleIO/HDFParticleIO.cpp 
)

SET(MESSAGE
  Message/Communicate.cpp 
)

SET(HAMILTONIAN
  QMCHamiltonians/QMCHamiltonian.cpp
)

#gsl library will be optional
IF(HAVE_LIBGSL)
  SET(UTILITIES ${UTILITIES} Optimize/GSLMinimize.cpp)
ENDIF(HAVE_LIBGSL)

SET(WFSBASE
  QMCWaveFunctions/JastrowBuilder.cpp
  QMCWaveFunctions/TrialWaveFunction.cpp
)

#Only basic packages common for both AtomicHF and QMC
SET(BASEDIR
   ${UTILITIES}
   ${MESSAGE}
)
ADD_LIBRARY(qmcbase ${BASEDIR})

#QMC-related package, build qmc
SET(QMCBASEDIR
   ${PARTICLE}
   ${PARTICLEIO}
   ${HAMILTONIAN}
   ${WFSBASE}
)

ADD_LIBRARY(qmc ${QMCBASEDIR})

IF(HAVE_OOMPI)
  CONFIGURE_FILE(${qmcPlusPlus_SOURCE_DIR}/src/oompi-config.h.cmake.in 
                 ${qmcPlusPlus_BINARY_DIR}/src/oompi-config.h)
  INCLUDE_DIRECTORIES(${qmcPlusPlus_SOURCE_DIR}/src/OOMPI)
  SUBDIRS(OOMPI)
  LINK_LIBRARIES(oompi)
ENDIF(HAVE_OOMPI)

SUBDIRS(QMCApp)

#do not build qmctools with mpi
IF(NOT HAVE_MPI)
  SUBDIRS(QMCTools)
ENDIF(NOT HAVE_MPI)
SUBDIRS(QMCTools)
