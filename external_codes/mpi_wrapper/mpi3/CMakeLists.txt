cmake_minimum_required(VERSION 3.16)

project(
  bmpi3
  VERSION 0.79.0
  DESCRIPTION
    "B-MPI3 is a C++ library wrapper for version 3.1 of the MPI standard interface that simplifies the utilization and maintenance of MPI code."
  HOMEPAGE_URL "https://gitlab.com/correaa/boost-mpi3"
  LANGUAGES CXX)

find_package(MPI REQUIRED COMPONENTS CXX) # C)  # might need to `module load
                                          # mpi`

add_library(${PROJECT_NAME} INTERFACE)

target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
            $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include/mpi3/dummy>)
target_link_libraries(${PROJECT_NAME} INTERFACE MPI::MPI_CXX)

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

# this makes CM FetchContent friendly
# https://www.foonathan.net/2022/06/cmake-fetchcontent/
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  return()
endif()

include(GNUInstallDirs)

include(CTest)

enable_testing()

add_subdirectory(include/mpi3/adaptors/fftw)

add_subdirectory(test)
