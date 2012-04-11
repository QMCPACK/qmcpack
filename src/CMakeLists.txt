CONFIGURE_FILE(${qmcpack_SOURCE_DIR}/src/config.h.cmake.in
               ${qmcpack_BINARY_DIR}/src/config.h)

####################################
# create libqmcutil
####################################
SET(UTILITIES
  Utilities/RandRandom.cpp 
  Utilities/RandomGenerator.cpp 
  Utilities/OhmmsObject.cpp 
  Utilities/SpeciesSet.cpp 
  Utilities/SimpleParser.cpp 
  Utilities/OhmmsInform.cpp 
  Utilities/OhmmsInfo.cpp 
  Utilities/NewTimer.cpp
  Utilities/ProgressReportEngine.cpp
  OhmmsData/Libxml2Doc.cpp
  OhmmsApp/ProjectData.cpp
  OhmmsApp/RandomNumberControl.cpp
  Numerics/OhmmsBlas.cpp
  Platforms/sysutil.cpp
  QMCFactory/OneDimGridFactory.cpp
  Message/Communicate.cpp 
  Message/MPIObjectBase.cpp 
  Optimize/VariableSet.cpp
  io/hdf_archive.cpp
)

IF(QMC_CUDA)
   SET(UTILITIES ${UTILITIES} 
       Numerics/CUDA/cuda_inverse.cu 
       CUDA/gpu_vector.cpp
       Particle/accept_kernel.cu)
ENDIF(QMC_CUDA)

IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(qmcutil ${UTILITIES})
ELSE(QMC_CUDA)
  ADD_LIBRARY(qmcutil ${UTILITIES})
ENDIF(QMC_CUDA)

#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcutil STATIC ${UTILITIES})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcutil SHARED ${UTILITIES})
#ENDIF(QMC_BUILD_STATIC)

####################################
# create liboompi
####################################
IF(HAVE_OOMPI)
  CONFIGURE_FILE(${qmcpack_SOURCE_DIR}/src/oompi-config.h.cmake.in 
                 ${qmcpack_BINARY_DIR}/src/oompi-config.h)
  INCLUDE_DIRECTORIES(${qmcpack_SOURCE_DIR}/src/OOMPI)
  SUBDIRS(OOMPI)

  set(MPI_LIBRARY oompi ${MPI_LIBRARY})

ENDIF(HAVE_OOMPI)

IF(OHMMS_DIM MATCHES 3)
  SET(PARTICLE Lattice/Uniform3DGridLayout.cpp)
ENDIF(OHMMS_DIM MATCHES 3)

IF(OHMMS_DIM MATCHES 2)
  SET(PARTICLE Lattice/Uniform2DGridLayout.cpp)
ENDIF(OHMMS_DIM MATCHES 2)

####################################
# create libqmcbase 
####################################
SET(PARTICLE ${PARTICLE}
  Particle/ParticleSet.cpp 
  Particle/ParticleSet.BC.cpp 
  Particle/MCWalkerConfiguration.cpp 
  Particle/DistanceTable.cpp
  Particle/HDFWalkerInputManager.cpp
  Particle/make_clones.cpp
  LongRange/KContainer.cpp
  LongRange/StructFact.cpp
  LongRange/LPQHIBasis.cpp
  LongRange/EwaldHandler.cpp
  LongRange/LRCoulombSingleton.cpp
)

IF(OHMMS_DIM MATCHES 2)
  SET(PARTICLE ${PARTICLE} LongRange/TwoDEwaldHandler.cpp)
ENDIF(OHMMS_DIM MATCHES 2)

IF(QMC_BUILD_LEVEL GREATER 2)
  IF(NOT QMC_COMPLEX)
    SET(PARTICLE ${PARTICLE} Particle/Bead_ParticleSet.cpp )
  ENDIF(NOT QMC_COMPLEX)
ENDIF(QMC_BUILD_LEVEL GREATER 2)

SET(PARTICLEIO
  ParticleTags.cpp 
  ParticleIO/ParticleLayoutIO.cpp 
  ParticleIO/XMLParticleIO.cpp 
  ParticleIO/ParticleIOUtility.cpp 
)

IF(HAVE_LIBHDF5)
  SET(PARTICLEIO ${PARTICLEIO}
      Particle/HDFWalkerIOEngine.cpp 
      Particle/HDFWalkerOutput.cpp
      Particle/HDFWalkerInput_0_0.cpp
      Particle/HDFWalkerInput_0_4.cpp
     )

   IF(OHMMS_DIM  MATCHES 3)
     SET(PARTICLEIO ${PARTICLEIO} ParticleIO/ESHDFParticleParser.cpp)
   ENDIF(OHMMS_DIM  MATCHES 3)

ENDIF(HAVE_LIBHDF5)

ADD_LIBRARY(qmcbase ${PARTICLE}  ${PARTICLEIO})

#IF(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcbase STATIC ${PARTICLE}  ${PARTICLEIO})
#ELSE(QMC_BUILD_STATIC)
#  ADD_LIBRARY(qmcbase SHARED ${PARTICLE}  ${PARTICLEIO})
#ENDIF(QMC_BUILD_STATIC)

SUBDIRS(QMCWaveFunctions)
SUBDIRS(QMCHamiltonians)
SUBDIRS(QMCDrivers)
SUBDIRS(QMCApp)

#SUBDIRS(einspline)


#IF(BUILD_QMCTOOLS)
  SUBDIRS(QMCTools)
#ENDIF(BUILD_QMCTOOLS)

##IF(BUILD_SANDBOX)
#  SUBDIRS(SandBox benchmark)
##ENDIF(BUILD_SANDBOX)

IF(BUILD_SQD)
  SUBDIRS(SQD)
ENDIF(BUILD_SQD)
