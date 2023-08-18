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
    MCCoords.cpp
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
    LongRange/EwaldHandlerQuasi2D.cpp
    LongRange/EwaldHandler3D.cpp
    LongRange/EwaldHandler2D.cpp
    LongRange/LRCoulombSingleton.cpp)

set(PARTICLEIO ParticleTags.cpp ParticleIO/LatticeIO.cpp ParticleIO/XMLParticleIO.cpp HDFWalkerOutput.cpp
               HDFWalkerInput_0_4.cpp)

if(USE_OBJECT_TARGET)
  add_library(qmcparticle OBJECT ${PARTICLE} ${PARTICLEIO})
else(USE_OBJECT_TARGET)
  add_library(qmcparticle ${PARTICLE} ${PARTICLEIO})
endif(USE_OBJECT_TARGET)

target_include_directories(qmcparticle PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcparticle PRIVATE platform_cpu_LA)
target_link_libraries(qmcparticle PUBLIC qmcnumerics qmcutil platform_runtime)
set(PARTICLE_OMPTARGET_SRCS
    createDistanceTableAAOMPTarget.cpp
    createDistanceTableABOMPTarget.cpp)

if(USE_OBJECT_TARGET)
  add_library(qmcparticle_omptarget OBJECT ${PARTICLE_OMPTARGET_SRCS})
else(USE_OBJECT_TARGET)
  add_library(qmcparticle_omptarget ${PARTICLE_OMPTARGET_SRCS})
endif(USE_OBJECT_TARGET)

target_include_directories(qmcparticle_omptarget PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcparticle_omptarget PUBLIC qmcutil containers)

target_compile_options(qmcparticle_omptarget PRIVATE "$<$<BOOL:${ENABLE_OFFLOAD_CLANG_DEBUG_O3}>:$<$<CONFIG:DEBUG>:-O3>>")

target_link_libraries(qmcparticle PUBLIC qmcparticle_omptarget)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
  add_subdirectory(Lattice/tests)
  add_subdirectory(LongRange/tests)
  add_subdirectory(ParticleIO/tests)
  add_subdirectory(ParticleBase/tests)
endif()
