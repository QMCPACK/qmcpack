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

SET(GITREV_BARE_FILE git-rev.h)
SET(GITREV_BARE_TMP git-rev-tmp.h)
SET(GITREV_FILE ${CMAKE_BINARY_DIR}/src/${GITREV_BARE_FILE})
SET(GITREV_TMP ${CMAKE_BINARY_DIR}/src/${GITREV_BARE_TMP})

# Always clean up git-rev-tmp.h file when cmake is run
EXECUTE_PROCESS(
  COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP}
)

IF (IS_GIT_PROJECT)
  # The following custom command picks up changes to the git revision information
  # every time the project is rebuilt. Even if the repositiory is updated (git pull)
  # without re-running cmake. It also appends '-dirty' to the commit hash if there are
  # unsaved changes to the repository.
  #
  # To avoid triggering a relink every time, the repository info is saved to
  # a temporary file, and the temporary file is copied over the actual file
  # only if the contents changed (using 'copy_if_different').
  # The temporary file is deleted to force the custom command to run on
  # the next build.
  #
  # Apparently custom commands need to be defined where the output is used.
  # If this in the main CMakeLists.txt it does not work.
  ADD_CUSTOM_COMMAND(
    OUTPUT ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_BRANCH_RAW " > ${GITREV_TMP}
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_HASH_RAW " >> ${GITREV_TMP}
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=40 --match="NoTagWithThisName" >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_COMMIT_LAST_CHANGED_RAW " >> ${GITREV_TMP}
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%ad >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo_append "#define GIT_COMMIT_SUBJECT_RAW \"" >> ${GITREV_TMP}
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%s | sed -E "s/\"/\\\\\"/g" | tr -d '\\n' >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo_append "\"" >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GITREV_TMP} ${GITREV_FILE}
    COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
  )

  # Print some configure-time git info (useful for understanding what commits
  # are in particular build for the nightly CDash reports)

  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE GIT_CONFIG_BRANCH
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Breaking down the arguments to 'git describe'
  #  --abbrev=40     Size of hash to print.  This should print the entire hash.
  #  --dirty         Append hash with '-dirty' if there are uncommited changes.
  # The behavior of describe looks for a tag in the parents first, and then falls
  # back to the commit hash (if --always is specified)
  #  --always        Show the commit hash as fallback
  #  --match="NoTagWithThisName"
  #     If a tag is found, the output looks like:
  #       second_annotated_tag-29-g1fd38cccc0fd2f683ec223ca0783bb671bfedd4e
  #     In order to always get just the commit hash, specify a tag pattern
  #     that should never match.
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=40 --match="NoTagWithThisName"
    OUTPUT_VARIABLE GIT_CONFIG_COMMIT_HASH
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  MESSAGE("Git branch: ${GIT_CONFIG_BRANCH}")
  MESSAGE("Git commit hash: ${GIT_CONFIG_COMMIT_HASH}")
ENDIF()


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
    Utilities/InfoStream.cpp
    Utilities/OutputManager.cpp
    Utilities/NewTimer.cpp
    Utilities/RunTimeManager.cpp
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
    spline2/einspline_allocator.c
    spline2/bspline_allocator.cpp
    spline2/MultiBspline.cpp
    )
IF (IS_GIT_PROJECT)
  SET(UTILITIES "${UTILITIES};${GITREV_TMP}")
ENDIF()

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
    Particle/VirtualParticleSet.cpp
    Particle/ParticleSet.BC.cpp 
    Particle/MCWalkerConfiguration.cpp 
    Particle/DistanceTableAA.cpp
    Particle/DistanceTableAB.cpp
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
      )
    IF(HAVE_ADIOS)
      SET(PARTICLEIO ${PARTICLEIO}
        ADIOS/ADIOS_verify.cpp
      )
    ENDIF(HAVE_ADIOS)

    IF(OHMMS_DIM  MATCHES 3)
      SET(PARTICLEIO ${PARTICLEIO} ParticleIO/ESHDFParticleParser.cpp)
    ENDIF(OHMMS_DIM  MATCHES 3)

  ENDIF(HAVE_LIBHDF5 OR HAVE_ADIOS)

  #Need to add this if nothing else but miniapps is compiled
  if(QMC_BUILD_LEVEL GREATER 4)
    set(PARTICLE  ${PARTICLE} QMCWaveFunctions/OrbitalBase.cpp)
  endif()

  ADD_LIBRARY(qmcbase ${PARTICLE}  ${PARTICLEIO}) 

  #IF(QMC_BUILD_STATIC)
  #  ADD_LIBRARY(qmcbase STATIC ${PARTICLE}  ${PARTICLEIO})
  #ELSE(QMC_BUILD_STATIC)
  #  ADD_LIBRARY(qmcbase SHARED ${PARTICLE}  ${PARTICLEIO})
  #ENDIF(QMC_BUILD_STATIC)
  #SUBDIRS(einspline)

  if(QMC_BUILD_LEVEL GREATER 4)

    MESSAGE(STATUS "Minimal build: only miniapps")
    SUBDIRS(miniapps)

  else() #{{{

  #############################
  # Creating ADIOS
  #############################

  SET(ADIOS ADIOS/ADIOS_config.cpp)
  IF(HAVE_ADIOS)
    SET(ADIOS ${ADIOS}
    Particle/AdiosWalkerInput.cpp
    ADIOS/ADIOS_profile.cpp
    ADIOS/ADIOS_verify.cpp
    ADIOS/ADIOSTrace.cpp
    )
  ENDIF(HAVE_ADIOS)

  ADD_LIBRARY(adios_config ${ADIOS})


  SUBDIRS(QMCWaveFunctions)
  SUBDIRS(QMCHamiltonians)
  SUBDIRS(QMCDrivers)

  IF(BUILD_AFQMC)
   SUBDIRS(AFQMC)
  ENDIF(BUILD_AFQMC)

  IF(BUILD_FCIQMC)
   SUBDIRS(FCIQMC)
  ENDIF(BUILD_FCIQMC)

  SUBDIRS(QMCApp)

  #IF(BUILD_QMCTOOLS)
  SUBDIRS(QMCTools)
  #ENDIF(BUILD_QMCTOOLS)

  IF(BUILD_SANDBOX)
    SUBDIRS(SandBox)
  ENDIF(BUILD_SANDBOX)

  IF(BUILD_SQD)
    SUBDIRS(SQD)
  ENDIF(BUILD_SQD)

  if (BUILD_UNIT_TESTS) #{
    SET(HAS_TARGET_COMPILE_DEFINITIONS 1)
    IF (CMAKE_VERSION VERSION_LESS "2.8.11")
      MESSAGE("CMake version 2.8.10 - some unit tests will not be built")
      SET(HAS_TARGET_COMPILE_DEFINITIONS 0)
    ENDIF()
    #Unit test directories
    INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/external_codes/catch)
    INCLUDE(${PROJECT_SOURCE_DIR}/CMake/unit_test.cmake)
    SUBDIRS(OhmmsPETE/tests)
    SUBDIRS(Numerics/tests)
    SUBDIRS(Utilities/tests)
    SUBDIRS(einspline/tests)
    SUBDIRS(OhmmsData/tests)
    SUBDIRS(io/tests)
    SUBDIRS(Lattice/tests)
    SUBDIRS(Particle/tests)
    SUBDIRS(ParticleIO/tests)
    SUBDIRS(QMCWaveFunctions/tests)
    SUBDIRS(QMCHamiltonians/tests)
    IF (HAS_TARGET_COMPILE_DEFINITIONS)
      SUBDIRS(ParticleBase/tests)
      SUBDIRS(Estimators/tests)
      SUBDIRS(QMCDrivers/tests)
      SUBDIRS(QMCApp/tests)
    ENDIF()
  endif() #}

endif() #}}}
