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
  LongRange/KContainer.cpp
  LongRange/StructFact.cpp
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


SET(OPTIMIZE
    Optimize/Minimize.cpp
    )

SET(WFSBASE
  QMCWaveFunctions/JastrowBuilder.cpp
  QMCWaveFunctions/TrialWaveFunction.cpp
)

#Only basic packages common for both AtomicHF and QMC
SET(BASEDIR
   ${UTILITIES}
   ${MESSAGE}
   ${OPTIMIZE}
)
ADD_LIBRARY(qmcbase ${BASEDIR})

#QMC-related package, build qmc
SET(QMCBASEDIR
   ${PARTICLE}
   ${PARTICLEIO}
   ${HAMILTONIAN}
   ${WFSBASE}
)

#SUBDIRS(SQD)
ADD_LIBRARY(qmc ${QMCBASEDIR})
SUBDIRS(QMCApp)
#SUBDIRS(QMC)
SUBDIRS(QMCTools)
