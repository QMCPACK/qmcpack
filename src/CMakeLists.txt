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

configure_file(${qmcpack_SOURCE_DIR}/src/config.h.cmake.in ${qmcpack_BINARY_DIR}/src/config.h)

include(${PROJECT_SOURCE_DIR}/CMake/unit_test.cmake)

set(GITREV_BARE_FILE git-rev.h)
set(GITREV_BARE_TMP git-rev-tmp.h)
set(GITREV_FILE ${qmcpack_BINARY_DIR}/src/${GITREV_BARE_FILE})
set(GITREV_TMP ${qmcpack_BINARY_DIR}/src/${GITREV_BARE_TMP})

# Always clean up git-rev-tmp.h file when cmake is run
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP})

if(IS_GIT_PROJECT)
  # The following custom command picks up changes to the git revision information
  # every time the project is rebuilt. Even if the repository is updated (git pull)
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

  # Sed flags were once an issue and some HPC have old sed's
  set(SED_FLAG "-E")
  execute_process(COMMAND "sed" ${SED_FLAG} "s/\"/\\\\\"/g" "<${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt" OUTPUT_QUIET
                  ERROR_VARIABLE SED_ERROR)
  if(SED_ERROR MATCHES ".*invalid.*")
    set(SED_FLAG "-r")
    execute_process(COMMAND "sed" ${SED_FLAG} "s/\"/\\\\\"/g" "<${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt" OUTPUT_QUIET
                            ERROR_QUIET ERROR_VARIABLE SED_ERROR)
    if(SED_ERROR MATCHES ".*invalid.*")
      message(
        WARNING
          "Your system supports neither the sed -E or -r flag, git revision information will not be included in output")
    else(SED_ERROR MATCHES ".*invalid.*")
      message("sed supports -r")
    endif(SED_ERROR MATCHES ".*invalid.*")
  else(SED_ERROR MATCHES ".*invalid.*")
    message("sed supports -E")
  endif(SED_ERROR MATCHES ".*invalid.*")

  add_custom_command(
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
    COMMAND ${GIT_EXECUTABLE} log -1 --format=%s | sed ${SED_FLAG} "s/\"/\\\\\"/g" | tr -d '\\n' >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo_append "\"" >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E echo >> ${GITREV_TMP}
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GITREV_TMP} ${GITREV_FILE}
    COMMAND ${CMAKE_COMMAND} -E remove ${GITREV_TMP}
    WORKING_DIRECTORY ${qmcpack_SOURCE_DIR}
    VERBATIM)

  # Print some configure-time git info (useful for understanding what commits
  # are in particular build for the nightly CDash reports)

  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    OUTPUT_VARIABLE GIT_CONFIG_BRANCH
    WORKING_DIRECTORY ${qmcpack_SOURCE_DIR} OUTPUT_STRIP_TRAILING_WHITESPACE)

  # Breaking down the arguments to 'git describe'
  #  --abbrev=40     Size of hash to print.  This should print the entire hash.
  #  --dirty         Append hash with '-dirty' if there are uncommitted changes.
  # The behavior of describe looks for a tag in the parents first, and then falls
  # back to the commit hash (if --always is specified)
  #  --always        Show the commit hash as fallback
  #  --match="NoTagWithThisName"
  #     If a tag is found, the output looks like:
  #       second_annotated_tag-29-g1fd38cccc0fd2f683ec223ca0783bb671bfedd4e
  #     In order to always get just the commit hash, specify a tag pattern
  #     that should never match.
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --always --dirty --abbrev=40 --match="NoTagWithThisName"
    OUTPUT_VARIABLE GIT_CONFIG_COMMIT_HASH
    WORKING_DIRECTORY ${qmcpack_SOURCE_DIR} OUTPUT_STRIP_TRAILING_WHITESPACE)
  message("Git branch: ${GIT_CONFIG_BRANCH}")
  message("Git commit hash: ${GIT_CONFIG_COMMIT_HASH}")
endif()

include_directories("${qmcpack_SOURCE_DIR}/external_codes/boost_multi")

set(RES_SRC
  Utiltites/ResourceCollection.cpp)

####################################
# create libqmcutil
####################################
set(UTILITIES
    Utilities/qmc_common.cpp
    Utilities/RandomGenerator.cpp
    Utilities/OhmmsObject.cpp
    Utilities/SpeciesSet.cpp
    Utilities/SimpleParser.cpp
    Utilities/Clock.cpp
    Utilities/NewTimer.cpp
    Utilities/TimerManager.cpp
    Utilities/RunTimeManager.cpp
    Utilities/ProgressReportEngine.cpp
    Utilities/unit_conversion.cpp
    Utilities/ResourceCollection.cpp
    OhmmsApp/ProjectData.cpp
    OhmmsApp/RandomNumberControl.cpp
    Numerics/OptimizableFunctorBase.cpp
    Numerics/SmoothFunctions.cpp
    Numerics/OneDimGridFactory.cpp
    Optimize/VariableSet.cpp)
if(IS_GIT_PROJECT)
  set(UTILITIES "${UTILITIES};${GITREV_TMP}")
endif()

if(USE_OBJECT_TARGET)
  add_library(qmcutil OBJECT ${UTILITIES})
else(USE_OBJECT_TARGET)
  add_library(qmcutil ${UTILITIES})
endif(USE_OBJECT_TARGET)
target_include_directories(qmcutil PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Utilities")
target_link_libraries(qmcutil PUBLIC message containers qmcio)
target_link_libraries(qmcutil PUBLIC LibXml2::LibXml2 Boost::boost ${QMC_UTIL_LIBS})

# For unit tests, enable use for the fake RNG
function(USE_FAKE_RNG TARGET)
  target_compile_definitions(${TARGET} PRIVATE "USE_FAKE_RNG")
endfunction()

# Put the fake RNG in a separate library so production code doesn't
# accidentally link to it
set(UNIT_TEST_UTILITIES Utilities/FakeRandom.cpp)
add_library(qmcfakerng ${UNIT_TEST_UTILITIES})

add_subdirectory(io)
add_subdirectory(einspline)
add_subdirectory(Containers)
add_subdirectory(Message)
add_subdirectory(Platforms)
add_subdirectory(Particle)
add_subdirectory(Utilities)

if(QMC_BUILD_SANDBOX_ONLY)

  message(STATUS "Minimal build: only Sandbox")
  add_subdirectory(Sandbox)

else() #{{{

  add_subdirectory(QMCWaveFunctions)
  add_subdirectory(QMCHamiltonians)
  add_subdirectory(Estimators)
  add_subdirectory(QMCDrivers)

  if(BUILD_LMYENGINE_INTERFACE)
    set(FORMIC_BINARY_DIR ${qmcpack_BINARY_DIR})
    add_subdirectory(formic/utils)
  endif(BUILD_LMYENGINE_INTERFACE)

  if(BUILD_AFQMC)
    add_subdirectory(AFQMC)
  endif(BUILD_AFQMC)

  if(BUILD_FCIQMC)
    add_subdirectory(FCIQMC)
  endif(BUILD_FCIQMC)

  add_subdirectory(QMCApp)
  add_subdirectory(QMCTools)

  if(BUILD_UNIT_TESTS) #{
    #Unit test directories
    add_subdirectory(Numerics/tests)
    add_subdirectory(Utilities/tests)
    add_subdirectory(einspline/tests)
    add_subdirectory(spline2/tests)
    add_subdirectory(Concurrency/tests)
    add_subdirectory(QMCWaveFunctions/tests)
    add_subdirectory(QMCHamiltonians/tests)
    add_subdirectory(type_traits/tests)
    add_subdirectory(QMCDrivers/tests)
    add_subdirectory(QMCApp/tests)
    add_subdirectory(Message/tests)
    add_subdirectory(Optimize/tests)
    add_subdirectory(QMCTools/tests)
  endif() #}

endif() #}}}
