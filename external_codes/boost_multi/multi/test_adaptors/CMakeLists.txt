# Copyright 2018-2024 Alfredo A. Correa
# Copyright 2024 Matt Borland
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt

cmake_minimum_required(VERSION 3.16)

# Library doesn't require installation, to still install this project:
# ~~~
# $ cmake .. --install-prefix=$HOME && cmake --build . --config Release --target test --target install -- -j $(nproc)
# ~~~
# to use this library in another CMake project
# ~~~
# project("Your project")
# find_package(boost-multi CONFIG REQUIRED)
# add_executable(${PROJECT_NAME} src/your_main.cpp)
# target_link_libraries(${PROJECT_NAME} boost-multi::boost-multi)
# ~~~

if (DEFINED BOOST_SUPERPROJECT_VERSION AND NOT DEFINED BOOST_MULTI_STANDALONE)

  project(boost_multi VERSION "${BOOST_SUPERPROJECT_VERSION}" LANGUAGES CXX)

  add_library(boost_multi INTERFACE)

  add_library(Boost::multi ALIAS boost_multi)

  target_include_directories(boost_multi INTERFACE include)

  target_compile_features(boost_multi INTERFACE cxx_std_14)

  message(STATUS "Boost.Multi: standalone mode OFF")

  if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")

    add_subdirectory(test)

  endif()

else()

  project(
    multi
    HOMEPAGE_URL "https://gitlab.com/correaa/boost-multi"
    DESCRIPTION "A header-only modern C++ library that provides access and manipulation of data in multidimensional arrays."
    VERSION 0.80.1
    LANGUAGES CXX)

  message(STATUS "Boost.Multi: standalone mode ON")

  message("current binary directory: ${CMAKE_CURRENT_BINARY_DIR}")

  include_directories(${PROJECT_SOURCE_DIR}/include)  # workaround for vscode to detect headers https://stackoverflow.com/a/68139743/225186

  include(CMakePackageConfigHelpers)
  include(CMakeDependentOption)
  include(GNUInstallDirs)

  find_program(MEMORYCHECK_COMMAND valgrind)
  set(VALGRIND_COMMAND_OPTIONS "-q --tool=memcheck --leak-check=yes --num-callers=51 --trace-children=yes --leak-check=full --track-origins=yes --gen-suppressions=all")

  # must go before `include(CTest)`
  set(MEMORYCHECK_COMMAND_OPTIONS "-q --tool=memcheck --leak-check=yes --num-callers=52 --trace-children=yes --leak-check=full --track-origins=yes --gen-suppressions=all")
  set(MEMORYCHECK_SUPPRESSIONS_FILE "${PROJECT_SOURCE_DIR}/.valgrind-suppressions")

  include(CTest)

  option(MULTI_BUILD_PACKAGE "Build package files as well" ON)

  cmake_dependent_option(MULTI_BUILD_TESTS "Enable multi tests" ON "BUILD_TESTING" OFF)

  cmake_dependent_option(MULTI_BUILD_PACKAGE_DEB "Create a DEB" ON "MULTI_BUILD_PACKAGE" OFF)

  add_library(multi INTERFACE)
  #target_compile_features(multi PUBLIC cxx_std_17)

  target_include_directories(multi INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
  target_compile_options(multi INTERFACE $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>: --expt-relaxed-constexpr --extended-lambda>)

  if(NOT CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    add_library(correaa::multi ALIAS multi)
  endif()

  # https://youtrack.jetbrains.com/issue/CPP-25608
  target_compile_features(${PROJECT_NAME} INTERFACE $<$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>:cxx_std_17>)
  target_compile_options(${PROJECT_NAME} INTERFACE $<$<CUDA_COMPILER_ID:NVIDIA>:-std=c++17>)

  # this makes CM FetchContent friendly https://www.foonathan.net/2022/06/cmake-fetchcontent/
  if(NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    return()
  endif()

  # Installation help
  configure_package_config_file("${PROJECT_SOURCE_DIR}/cmake/multi-config.cmake.in" "${PROJECT_BINARY_DIR}/multi-config.cmake" INSTALL_DESTINATION "share/cmake/multi")

  write_basic_package_version_file("${PROJECT_BINARY_DIR}/multi-config-version.cmake" COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)

  message("current install prefix directory: ${CMAKE_INSTALL_PREFIX}")

  install(
    TARGETS multi
    EXPORT multi-targets
    INCLUDES
    DESTINATION "${CMAKE_INSTALL_DATADIR}")

  install(
    EXPORT ${PROJECT_NAME}-targets
    DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}"
    NAMESPACE multi::
    FILE "${PROJECT_NAME}-targets.cmake")

  install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake" DESTINATION "${CMAKE_INSTALL_DATADIR}/cmake/${PROJECT_NAME}")

  install(DIRECTORY "include/" TYPE INCLUDE)

  add_subdirectory(test)

  # if(MULTI_BUILD_TESTS)
    find_package(Boost 1.65 COMPONENTS)  # unit_test_framework)  # 1.65 needed for BOOST_TEST_GLOBAL_FIXTURE, you can use your own Boost and use `cmake -DBOOST_ROOT=$HOME/local`
    if(NOT Boost_FOUND)
      message(WARNING "Cannot find Boost 1.65+, Multi library will not be fully tested.")
    else()
      enable_testing()

      add_subdirectory(include/boost/multi/adaptors/blas)
      add_subdirectory(include/boost/multi/adaptors/complex)
      add_subdirectory(include/boost/multi/adaptors/cuda)
      add_subdirectory(include/boost/multi/adaptors/fftw)

      find_package(LAPACK)
      if(LAPACK_FOUND)
        add_subdirectory(include/boost/multi/adaptors/lapack)
      endif()

      add_subdirectory(include/boost/multi/adaptors/thrust)
      if(ENABLE_CUDA)
        add_subdirectory(include/boost/multi/adaptors/cufft)
      endif()
      if(ENABLE_HIP)
        add_subdirectory(include/boost/multi/adaptors/hipfft)
        add_subdirectory(include/boost/multi/adaptors/hipthrust/test)
      endif()
    endif()
  # endif()

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

endif()
