cmake_minimum_required(VERSION 3.16)

project(
  bmpi3
  VERSION 0.79.0
  # DESCRIPTION
  #   "B-MPI3 is a C++ library wrapper for version 3.1 of the MPI standard interface that simplifies the utilization and maintenance of MPI code."
  # HOMEPAGE_URL "https://gitlab.com/correaa/boost-mpi3"
  LANGUAGES CXX)

include(GNUInstallDirs)

# Find the dependencies.
find_package(MPI REQUIRED COMPONENTS CXX) # C)  # might need to `module load
                                          # mpi`
find_package(Boost REQUIRED COMPONENTS serialization)

# Set up the (interface) library.
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/mpi3/dummy>)
target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_CXX)
target_link_libraries(${PROJECT_NAME} INTERFACE Boost::serialization)

# Set up the installation of the library/project.
option(BMPI3_INSTALL "Install the bmpi3 library with the current project" ON)
if(BMPI3_INSTALL)
  include(CMakePackageConfigHelpers)
  install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets)
  install(
    DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  install(
    EXPORT ${PROJECT_NAME}Targets
    FILE "${PROJECT_NAME}-config-targets.cmake"
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
  write_basic_package_version_file( "${PROJECT_NAME}-config-version.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY AnyNewerVersion)
  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
endif()

# ~~~
# to use this project directly from CMake
# FetchContent_Declare(
#     bmpi3
#     GIT_REPOSITORY         git@gitlab.com:correaa/boost-mpi3.git  # https://gitlab.com/correaa/boost-mpi3.git
#     GIT_TAG                master)
# FetchContent_MakeAvailable(bmpi3)
# add_executable(main main.cpp)
# target_link_libraries(main PUBLIC bmpi3)
# ~~~

# Set up the tests.
set(BMPI3_BUILD_TESTING_DEFAULT OFF)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(BMPI3_BUILD_TESTING_DEFAULT ON)
endif()
option(BUILD_TESTING "Enable testing with CTest."
  ${BMPI3_BUILD_TESTING_DEFAULT})
option(BMPI3_BUILD_TESTING "Enable testing bmpi3 with CTest."
  ${BMPI3_BUILD_TESTING_DEFAULT})
include(CTest)

if(BUILD_TESTING AND BMPI3_BUILD_TESTING)
  find_package(PkgConfig)
  add_subdirectory(include/mpi3/adaptors/fftw)
  add_subdirectory(test)
endif()
