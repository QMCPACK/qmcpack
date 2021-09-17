#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(UTILITIES
    qmc_common.cpp
    RandomGenerator.cpp
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
target_link_libraries(qmcutil PUBLIC message qmcio)
target_link_libraries(qmcutil PUBLIC LibXml2::LibXml2 Boost::boost ${QMC_UTIL_LIBS})

# Put the fake RNG in a separate library so production code doesn't
# accidentally link to it
set(UNIT_TEST_UTILITIES FakeRandom.cpp)
add_library(qmcfakerng ${UNIT_TEST_UTILITIES})

subdirs(for_testing)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
