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
  Platforms/sysutil.cpp
)

SET(PARTICLE
  Particle/ParticleSet.cpp 
  Particle/MCWalkerConfiguration.cpp 
  Particle/DistanceTable.cpp
  Particle/HDFWalkerIO.cpp
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

SET(OPTIMIZE
    Optimize/Minimize.cpp
    )

SET(HAMILTONIAN
  QMCHamiltonians/QMCHamiltonianBase.cpp 
  QMCHamiltonians/LocalPPotential.cpp
)

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
   ${OPTIMIZE}
   ${HAMILTONIAN}
   ${WFSBASE}
)

SUBDIRS(SQD PseudoGen)
#ADD_LIBRARY(qmc ${QMCBASEDIR})
#SUBDIRS(QMC)

