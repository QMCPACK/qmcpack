#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////


####################################
# create libqmcparticle
####################################
SET(PARTICLE
  InitMolecularSystem.cpp
  ParticleSetPool.cpp
  ParticleSet.cpp
  PSdispatcher.cpp
  VirtualParticleSet.cpp
  ParticleSet.BC.cpp
  DynamicCoordinatesBuilder.cpp
  MCWalkerConfiguration.cpp
  WalkerConfigurations.cpp
  SampleStack.cpp
  createDistanceTableAA.cpp
  createDistanceTableAB.cpp
  HDFWalkerInputManager.cpp
  LongRange/KContainer.cpp
  LongRange/StructFact.cpp
  LongRange/LPQHIBasis.cpp
  LongRange/LPQHISRCoulombBasis.cpp
  LongRange/EwaldHandler.cpp
  LongRange/EwaldHandler3D.cpp
  LongRange/LRCoulombSingleton.cpp
  )

IF(ENABLE_OFFLOAD)
  SET(PARTICLE ${PARTICLE}
    createDistanceTableAAOMPTarget.cpp
    createDistanceTableABOMPTarget.cpp
    )
ENDIF(ENABLE_OFFLOAD)

IF(OHMMS_DIM MATCHES 2)
  SET(PARTICLE ${PARTICLE} LongRange/TwoDEwaldHandler.cpp)
ENDIF(OHMMS_DIM MATCHES 2)

SET(PARTICLEIO
  ParticleTags.cpp
  ParticleIO/ParticleLayoutIO.cpp
  ParticleIO/XMLParticleIO.cpp
  ParticleIO/ParticleIOUtility.cpp
  HDFWalkerOutput.cpp
  HDFWalkerInput_0_0.cpp
  HDFWalkerInput_0_4.cpp
  )

IF(OHMMS_DIM MATCHES 3)
  SET(PARTICLEIO ${PARTICLEIO} ParticleIO/ESHDFParticleParser.cpp)
ENDIF(OHMMS_DIM MATCHES 3)

IF(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmcparticle OBJECT ${PARTICLE} ${PARTICLEIO})
ELSE(USE_OBJECT_TARGET)
  ADD_LIBRARY(qmcparticle ${PARTICLE} ${PARTICLEIO})
ENDIF(USE_OBJECT_TARGET)

TARGET_INCLUDE_DIRECTORIES(qmcparticle PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
TARGET_LINK_LIBRARIES(qmcparticle PRIVATE Math::BLAS_LAPACK Math::scalar_vector_functions)
TARGET_LINK_LIBRARIES(qmcparticle PUBLIC qmcutil platform_device)

IF(QMC_CUDA)
  CUDA_ADD_LIBRARY(qmcparticle_cuda accept_kernel.cu)
  TARGET_LINK_LIBRARIES(qmcparticle_cuda PRIVATE platform_device)
  TARGET_LINK_LIBRARIES(qmcparticle PRIVATE qmcparticle_cuda)
ENDIF(QMC_CUDA)

IF(BUILD_UNIT_TESTS)
  SUBDIRS(tests)
  SUBDIRS(Lattice/tests)
  SUBDIRS(LongRange/tests)
  SUBDIRS(ParticleIO/tests)
  SUBDIRS(ParticleBase/tests)
ENDIF()
