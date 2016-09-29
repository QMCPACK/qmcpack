#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign  
#//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
#//                    Norbert Podhorszki, pnorbert@ornl.gov, Oak Ridge National Laboratory
#//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
#//                    Ying Wai Li, yingwaili@ornl.gov, Oak Ridge National Laboratory
#//                    Marl Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
#//		       Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////

CONFIGURE_FILE(${qmcpack_SOURCE_DIR}/src/config.h.cmake.in
  ${qmcpack_BINARY_DIR}/src/config.h)

if(QMC_PHI)
  MESSAGE(STATUS "Minimal build to test codelets on phi")
else()
  ####################################
  # create libqmcutil
  ####################################
  SET(UTILITIES
    qmc_common.cpp
    Utilities/RandRandom.cpp 
    Utilities/RandomGenerator.cpp 
    Utilities/OhmmsObject.cpp 
    Utilities/SpeciesSet.cpp 
    Utilities/SimpleParser.cpp 
    Utilities/OhmmsInform.cpp 
    Utilities/OhmmsInfo.cpp 
    Utilities/NewTimer.cpp
    Utilities/ProgressReportEngine.cpp
    Utilities/unit_conversion.cpp
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

  IF(QMC_ADIOS)
    SET(UTILITIES ${UTILITIES}
      ADIOS/ADIOS_verify.cpp
      ADIOS/ADIOS_config.h
      ADIOS/ADIOS_config.cpp
      )
  ENDIF(QMC_ADIOS)

  IF(QMC_CUDA)
    SET(UTILITIES ${UTILITIES} 
      Numerics/CUDA/cuda_inverse.cu 
      CUDA/gpu_vector.cpp
      CUDA/gpu_misc.cpp
      Particle/accept_kernel.cu)
  ENDIF(QMC_CUDA)

  IF(QMC_CUDA)
    CUDA_ADD_LIBRARY(qmcutil ${UTILITIES})
    CUDA_ADD_CUBLAS_TO_TARGET(qmcutil)
  ELSE(QMC_CUDA)
    ADD_LIBRARY(qmcutil ${UTILITIES} ${ADIOS})
  ENDIF(QMC_CUDA)
  
  # For unit tests, enable use for the fake RNG
  FUNCTION( USE_FAKE_RNG TARGET )
    TARGET_COMPILE_DEFINITIONS(${TARGET} PRIVATE "USE_FAKE_RNG")
  ENDFUNCTION()

  # Put the fake RNG in a separate library so production code doesn't
  # accidentally link to it
  SET(UNIT_TEST_UTILITIES
    Utilities/FakeRandom.cpp 
  ) 
  ADD_LIBRARY(qmcfakerng ${UNIT_TEST_UTILITIES})

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
    IF ( MPI_LIBRARY )
       SET(MPI_LIBRARY oompi ${MPI_LIBRARY})
    ELSE()
       SET(MPI_LIBRARY oompi )
    ENDIF()
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
    Particle/VirtualParticleSet.cpp 
    Particle/MCWalkerConfiguration.cpp 
    Particle/DistanceTable.cpp
    Particle/HDFWalkerInputManager.cpp
    Particle/make_clones.cpp
    LongRange/KContainer.cpp
    LongRange/StructFact.cpp
    LongRange/LPQHIBasis.cpp
    LongRange/LPQHISRCoulombBasis.cpp
    LongRange/EwaldHandler.cpp
    LongRange/LRCoulombSingleton.cpp
    )

  IF(OHMMS_DIM MATCHES 2)
    SET(PARTICLE ${PARTICLE} LongRange/TwoDEwaldHandler.cpp)
  ENDIF(OHMMS_DIM MATCHES 2)

  # IF(QMC_BUILD_LEVEL GREATER 2)
  #   IF(NOT QMC_COMPLEX)
  #     SET(PARTICLE ${PARTICLE} Particle/Bead_ParticleSet.cpp )
  #   ENDIF(NOT QMC_COMPLEX)
  # ENDIF(QMC_BUILD_LEVEL GREATER 2)

  SET(PARTICLEIO
    ParticleTags.cpp 
    ParticleIO/ParticleLayoutIO.cpp 
    ParticleIO/XMLParticleIO.cpp 
    ParticleIO/ParticleIOUtility.cpp 
    )

  IF(HAVE_LIBHDF5 OR HAVE_ADIOS)
    SET(PARTICLEIO ${PARTICLEIO}
      Particle/HDFWalkerIOEngine.cpp 
      Particle/HDFWalkerOutput.cpp
      Particle/HDFWalkerInput_0_0.cpp
      Particle/HDFWalkerInput_0_4.cpp
      ADIOS/ADIOS_verify.cpp
      )

    IF(OHMMS_DIM  MATCHES 3)
      SET(PARTICLEIO ${PARTICLEIO} ParticleIO/ESHDFParticleParser.cpp)
    ENDIF(OHMMS_DIM  MATCHES 3)

  ENDIF(HAVE_LIBHDF5 OR HAVE_ADIOS)

  ADD_LIBRARY(qmcbase ${PARTICLE}  ${PARTICLEIO}) 

  #IF(QMC_BUILD_STATIC)
  #  ADD_LIBRARY(qmcbase STATIC ${PARTICLE}  ${PARTICLEIO})
  #ELSE(QMC_BUILD_STATIC)
  #  ADD_LIBRARY(qmcbase SHARED ${PARTICLE}  ${PARTICLEIO})
  #ENDIF(QMC_BUILD_STATIC)

  #SUBDIRS(einspline)

  #############################
  # Creating ADIOS
  #############################


  SET(ADIOS ADIOS/ADIOS_config.cpp)
  #if(HAVE_ADIOS)
    SET(ADIOS ${ADIOS}
    Particle/AdiosWalkerInput.cpp
    ADIOS/ADIOS_profile.cpp
    ADIOS/ADIOS_verify.cpp
    ADIOS/ADIOSTrace.cpp
    )
  #endif(HAVE_ADIOS)

  ADD_LIBRARY(adios_config ${ADIOS})

  SUBDIRS(QMCWaveFunctions)
  SUBDIRS(QMCHamiltonians)
  SUBDIRS(QMCDrivers)
  SUBDIRS(QMCApp)

  #IF(BUILD_QMCTOOLS)
  SUBDIRS(QMCTools)
  #ENDIF(BUILD_QMCTOOLS)

  IF(BUILD_SANDBOX)
    SUBDIRS(SandBox)
  ENDIF(BUILD_SANDBOX)

  IF(BUILD_BENCHMARK)
    SUBDIRS(benchmark)
  ENDIF(BUILD_BENCHMARK)

  IF(BUILD_SQD)
    SUBDIRS(SQD)
  ENDIF(BUILD_SQD)

  if (BUILD_UNIT_TESTS)
    SET(HAS_TARGET_COMPILE_DEFINITIONS 1)
    IF (CMAKE_VERSION VERSION_LESS "2.8.11")
      MESSAGE("CMake version 2.8.10 - some unit tests will not be built")
      SET(HAS_TARGET_COMPILE_DEFINITIONS 0)
    ENDIF()
    #Unit test directories
    INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external_codes/catch)
    INCLUDE(${PROJECT_SOURCE_DIR}/CMake/unit_test.cmake)
    SUBDIRS(Numerics/tests)
    SUBDIRS(Utilities/tests)
    SUBDIRS(einspline/tests)
    SUBDIRS(OhmmsData/tests)
    SUBDIRS(Lattice/tests)
    SUBDIRS(Particle/tests)
    SUBDIRS(ParticleIO/tests)
    SUBDIRS(QMCWaveFunctions/tests)
    SUBDIRS(QMCHamiltonians/tests)
    IF (HAS_TARGET_COMPILE_DEFINITIONS)
      SUBDIRS(ParticleBase/tests)
      SUBDIRS(QMCDrivers/tests)
    ENDIF()
  endif()
endif()
