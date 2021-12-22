#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Peter Doak, , doakpw@ornl.gov, Oak Ridge National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

# there is a provides a utilities target for qmcio that doesn't create a circular
# dependence
add_library(qmcutil_basic ModernStringUtils.cpp)
target_include_directories(qmcutil_basic PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

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
add_library(qmcutil_integrated ${UTILITIES})

if(IS_GIT_PROJECT)
  add_dependencies(qmcutil_integrated gitrev)
endif()

target_link_libraries(qmcutil_integrated PUBLIC message qmcio)
target_link_libraries(qmcutil_integrated PUBLIC LibXml2::LibXml2 Boost::boost ${QMC_UTIL_LIBS})

add_library(qmcutil INTERFACE)
target_include_directories(qmcutil INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcutil INTERFACE qmcutil_basic qmcutil_integrated)


# Put the fake RNG in a separate library so production code doesn't
# accidentally link to it
set(UNIT_TEST_UTILITIES FakeRandom.cpp)
add_library(qmcfakerng ${UNIT_TEST_UTILITIES})

subdirs(for_testing)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
