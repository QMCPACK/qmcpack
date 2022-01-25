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
set(PARTICLE
    InitMolecularSystem.cpp
    SimulationCell.cpp
    ParticleSetPool.cpp
    ParticleSet.cpp
    PSdispatcher.cpp
    VirtualParticleSet.cpp
    ParticleSet.BC.cpp
    DynamicCoordinatesBuilder.cpp
    MCWalkerConfiguration.cpp
    WalkerConfigurations.cpp
    SpeciesSet.cpp
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
    LongRange/LRCoulombSingleton.cpp)

if(ENABLE_OFFLOAD)
  set(PARTICLE ${PARTICLE} createDistanceTableAAOMPTarget.cpp createDistanceTableABOMPTarget.cpp)
endif(ENABLE_OFFLOAD)

if(OHMMS_DIM MATCHES 2)
  set(PARTICLE ${PARTICLE} LongRange/TwoDEwaldHandler.cpp)
endif(OHMMS_DIM MATCHES 2)

set(PARTICLEIO
    ParticleTags.cpp
    ParticleIO/ParticleLayoutIO.cpp
    ParticleIO/XMLParticleIO.cpp
    HDFWalkerOutput.cpp
    HDFWalkerInput_0_4.cpp)

if(USE_OBJECT_TARGET)
  add_library(qmcparticle OBJECT ${PARTICLE} ${PARTICLEIO})
else(USE_OBJECT_TARGET)
  add_library(qmcparticle ${PARTICLE} ${PARTICLEIO})
endif(USE_OBJECT_TARGET)

target_include_directories(qmcparticle PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcparticle PRIVATE platform_cpu_LA)
target_link_libraries(qmcparticle PUBLIC qmcnumerics qmcutil platform_runtime)

if(QMC_CUDA)
  if(NOT QMC_CUDA2HIP)
    add_library(qmcparticle_cuda accept_kernel.cu)
  else()
    hip_add_library(qmcparticle_cuda accept_kernel.cu)
  endif(NOT QMC_CUDA2HIP)
  target_link_libraries(qmcparticle_cuda PRIVATE platform_runtime)
  target_link_libraries(qmcparticle PRIVATE qmcparticle_cuda)
endif(QMC_CUDA)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
  add_subdirectory(Lattice/tests)
  add_subdirectory(LongRange/tests)
  add_subdirectory(ParticleIO/tests)
  add_subdirectory(ParticleBase/tests)
endif()
