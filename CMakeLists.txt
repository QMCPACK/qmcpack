# Copyright 2018-2025 Alfredo A. Correa
# Copyright 2024 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.16)

if (DEFINED BOOST_SUPERPROJECT_VERSION AND NOT DEFINED BOOST_MULTI_STANDALONE)

  project(boost_multi VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX)

  add_library(boost_multi INTERFACE)

  add_library(Boost::multi ALIAS boost_multi)

  target_include_directories(boost_multi INTERFACE include)

  target_compile_features(boost_multi INTERFACE cxx_std_17)

  message(STATUS "Boost.Multi: standalone mode OFF")

  if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")

    add_subdirectory(test)

  endif()

  return()
endif()

# Library doesn't require installation, to still install this project:
# ~~~
# $ cmake .. --install-prefix=$HOME && cmake --build . --config Release --target test --target install -- --parallel
# ~~~
# to use this library in another CMake project
# ~~~
# project("Your project")
# find_package(boost-multi CONFIG REQUIRED)
# add_executable(${PROJECT_NAME} src/your_main.cpp)
# target_link_libraries(${PROJECT_NAME} boost-multi::boost-multi)
# ~~~

message("current binary directory: ${CMAKE_CURRENT_BINARY_DIR}")

project(
  multi
  HOMEPAGE_URL "https://gitlab.com/correaa/boost-multi"
  DESCRIPTION "A header-only modern C++ library that provides access and manipulation of data in multidimensional arrays."
  VERSION 0.88.0
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Boost.Multi: standalone mode ON")

add_library(multi INTERFACE)

add_library(boost.multi ALIAS multi)
add_library(boost::multi ALIAS multi)

#target_sources(
#  multi
#  PUBLIC
#    FILE_SET HEADERS
#    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
#    FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/boost/multi/array.hpp
#)

set_target_properties(
  multi
  PROPERTIES
    # VERIFY_INTERFACE_HEADER_SET ON
    EXPORT_NAME multi-targets
)

target_include_directories(multi INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

# Documentation driver for MrDocs — configured but not built by default.
# Produces a compile_commands.json entry for the public headers when
# cmake is invoked with -DCMAKE_EXPORT_COMPILE_COMMANDS=ON.
add_library(mrdocs_driver OBJECT ${PROJECT_SOURCE_DIR}/doc/mrdocs_driver.cpp)
target_include_directories(mrdocs_driver PRIVATE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)
target_compile_features(mrdocs_driver PRIVATE cxx_std_17)
set_target_properties(mrdocs_driver PROPERTIES EXCLUDE_FROM_ALL TRUE)

#if(NOT ENABLE_CUDA)
#target_compile_options(multi INTERFACE
#  $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>: --expt-relaxed-constexpr --extended-lambda>
#  $<$<CXX_COMPILER_ID:MSVC>:
#    /EHsc
#    $<$<COMPILE_FEATURES:cxx_std_17>:/permissive- >
#    /wd5027  # move assignment operator was implicitly defined as deleted, for detail::array_allocator
#  >  # this is needed for MSVC in C++17 mode
#)
#else()
#target_compile_options(multi INTERFACE
#  $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>: --expt-relaxed-constexpr --extended-lambda>
#  $<$<CXX_COMPILER_ID:MSVC>:-Xcompiler=/EHsc $<$<COMPILE_FEATURES:cxx_std_17>:-Xcompiler=/permissive- > >  # this is needed for MSVC in C++17 mode
#)

if(NOT ENABLE_CIRCLE)
  target_compile_options(multi INTERFACE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CUDA_COMPILER_ID:Clang>>:
      -Wno-c++17-extensions
      -Wno-c++98-compat -Wno-c++98-compat-bind-to-temporary-copy
      -Wno-padded
      $<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,21>:-Wno-nrvo>
      $<$<VERSION_GREATER_EQUAL:$<CXX_COMPILER_VERSION>,12>:-Wno-pre-c++17-compat>
    >
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>: --expt-relaxed-constexpr>  # --extended-lambda>
    $<$<CXX_COMPILER_ID:MSVC>:
      $<$<NOT:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:
        /EHsc /permissive-
      >
      $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
        -Xcompiler=/EHsc -Xcompiler=/permissive- -Xcompiler=/Zc:preprocessor
      >
    >
  )
endif()

# https://youtrack.jetbrains.com/issue/CPP-25608
target_compile_features(multi INTERFACE cxx_std_17)
if(ENABLE_CUDA)
  if(WIN32)
    target_compile_features(multi INTERFACE cuda_std_20)
  else()
    target_compile_features(multi INTERFACE cuda_std_17)
  endif()
endif()

if(NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  add_library(correaa::multi ALIAS multi)
endif()

# this makes CM FetchContent friendly https://www.foonathan.net/2022/06/cmake-fetchcontent/
if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  return()
endif()
################################################################################


include(CMakePackageConfigHelpers)
include(CMakeDependentOption)
include(GNUInstallDirs)

include(CTest)

option(MULTI_BUILD_PACKAGE "Build package files as well" ON)

cmake_dependent_option(MULTI_BUILD_TESTS "Enable multi tests" ON "BUILD_TESTING" OFF)

cmake_dependent_option(MULTI_BUILD_PACKAGE_DEB "Create a DEB" ON "MULTI_BUILD_PACKAGE" OFF)

# Installation help
configure_package_config_file("${PROJECT_SOURCE_DIR}/cmake/multi-config.cmake.in" "${PROJECT_BINARY_DIR}/multi-config.cmake" INSTALL_DESTINATION "share/cmake/multi")

write_basic_package_version_file("${PROJECT_BINARY_DIR}/multi-config-version.cmake" COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)

message("current install prefix directory: ${CMAKE_INSTALL_PREFIX}")

install(
  TARGETS multi COMPONENT multi
  EXPORT multi
  # EXPORT multi-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}$<$<CONFIG:Debug>:/debug>
  # "${CMAKE_INSTALL_DATADIR}"
  # RUNTIME_DESTINATION ${CMAKE_INSTALL_BINDIR}$<$<CONFIG:Debug>:/debug>
  # INCLUDES
  # FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

#install(
#  EXPORT ${PROJECT_NAME}-targets
#  DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}"
#  NAMESPACE multi::
#  FILE "${PROJECT_NAME}-targets.cmake")

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake" DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}")

install(DIRECTORY "include/" TYPE INCLUDE)

if(MULTI_BUILD_TESTS)
#  find_package(Boost QUIET COMPONENTS headers)
# find_package(Boost COMPONENTS NO_MODULE)  # you can use your own Boost and use `cmake -DBOOST_ROOT=$HOME/local`
  find_package(Boost CONFIG QUIET)  # you can use your own Boost with `cmake -DBoost_ROOT=$HOME/local`
  if(NOT Boost_FOUND)
    message(WARNING "Cannot find Boost, Multi library will not be fully tested.\n  Install with: sudo apt install --no-install-recommends libboost-dev")
  else()
    enable_testing()

    add_subdirectory(test)

    add_subdirectory(include/boost/multi/adaptors/blas)
    # # add_subdirectory(include/boost/multi/adaptors/complex)
    add_subdirectory(include/boost/multi/adaptors/fftw)

    find_package(LAPACK)
    if(LAPACK_FOUND)
      add_subdirectory(include/boost/multi/adaptors/lapack)
    endif()

    add_subdirectory(include/boost/multi/adaptors/thrust)
    if(ENABLE_CUDA)
      add_subdirectory(include/boost/multi/adaptors/cuda)
      add_subdirectory(include/boost/multi/adaptors/cufft)
      add_subdirectory(include/boost/multi/adaptors/cutensor)
    endif()
    if(ENABLE_HIP)
    #  add_subdirectory(include/boost/multi/adaptors/hipfft)
    #  add_subdirectory(include/boost/multi/adaptors/hipthrust/test)
    endif()

    add_subdirectory(include/boost/multi/adaptors/tblis)
    add_subdirectory(include/boost/multi/adaptors)
  endif()
endif()

# # Find the quom tool. Adjust the path if quom is not in your system's PATH.
# # Alternatively, you could download it as part of a custom build step.
# find_program(QUOM_EXECUTABLE quom REQUIRED)

# # Define variables for the input and output files
# set(QUOM_INPUT_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/boost/multi/array.hpp")
# set(QUOM_OUTPUT_HEADER "${CMAKE_CURRENT_BINARY_DIR}/array_amalgamated.hpp")

# # Use add_custom_command to run quom and generate the single header.
# add_custom_command(
#     OUTPUT ${QUOM_OUTPUT_HEADER}
#     COMMAND ${QUOM_EXECUTABLE} --output ${QUOM_OUTPUT_HEADER} ${QUOM_INPUT_HEADER}
#     DEPENDS ${QUOM_INPUT_HEADER}
#     COMMENT "Generating single header file with quom"
# )

# # Create an interface library to encapsulate the single header.
# # This makes it easy for other targets to use the generated header.
# add_library(my_project_single_header INTERFACE)

# # The generated header becomes a dependency of the interface library.
# add_dependencies(my_project_single_header ${QUOM_OUTPUT_HEADER})

# # Add the binary directory to the include path so targets can find the generated header.
# target_include_directories(my_project_single_header INTERFACE
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
# )

if(MULTI_BUILD_PACKAGE)
  list(APPEND source-generators TBZ2 TGZ TXZ ZIP)

  if(CMAKE_HOST_WIN32)
    list(APPEND binary-generators "WIX")
  endif()

  if(MULTI_BUILD_PACKAGE_DEB)
    list(APPEND binary-generators "DEB")
  endif()

  if(MULTI_BUILD_RPM)
    list(APPEND binary-generators "RPM")
  endif()

  set(CPACK_PACKAGE_NAME ${PROJECT_NAME}
    CACHE STRING "The resulting package name"
  )
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "all")
  set(CPACK_PACKAGE_ARCHITECTURE "all")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Multidimensional arrays for C++"
    CACHE STRING "Package description for the package metadata"
  )

  set(CPACK_PACKAGE_VENDOR "alfredo.correa@gmail.com")
  # set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})
  # SET(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_SOURCE_DIR}/_packages")
  # set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/some")#/${CMAKE_PROJECT_VERSION}")

  # set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
  # set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
  # set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})

  set(CPACK_PACKAGE_CONTACT "alfredo.correa@gmail.com")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Alfredo A. Correa")

  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
  set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")

  set(CPACK_DEBIAN_FILE_NAME "multi_all.deb")
  # set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE)#ONE_PER_GROUP)
  # set(CPACK_DEB_COMPONENT_INSTALL YES)

  set(CPACK_SOURCE_GENERATOR ${source-generators})
  set(CPACK_GENERATOR ${binary-generators})

  # set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
  # set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}")

  list(
    APPEND
    CPACK_SOURCE_IGNORE_FILES
    /.git/
    /.build*/
    /build/
    .gitignore
    .DS_Store)

  include(CPack)
endif()
