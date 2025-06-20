#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

add_library(cxx_helpers ModernStringUtils.cpp)
target_include_directories(cxx_helpers PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

set(QMC_RNG FakeRandom.cpp RandomGenerator.cpp StdRandom.cpp)
add_library(qmcrng ${QMC_RNG})

set(UTILITIES
    qmc_common.cpp
    OhmmsObject.cpp
    SimpleParser.cpp
    Clock.cpp
    NewTimer.cpp
    TimerManager.cpp
    RunTimeManager.cpp
    ProgressReportEngine.cpp
    unit_conversion.cpp
    ResourceCollection.cpp
    ProjectData.cpp
    RandomNumberControl.cpp)
add_library(qmcutil ${UTILITIES})

if(IS_GIT_PROJECT)
  add_dependencies(qmcutil gitrev)
endif()

target_include_directories(qmcutil PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcutil PUBLIC cxx_helpers message qmcio qmcrng)
target_link_libraries(qmcutil PUBLIC LibXml2::LibXml2 ${QMC_UTIL_LIBS})

subdirs(for_testing)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
