# -*-indent-tabs-mode:nil;c-basic-offset:2;tab-width:4;autowrap:nil;-*-
cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(boost-multi-adaptors-blas VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if((NOT CMAKE_CXX_COMPILER_ID STREQUAL "PGI") AND (NOT DART_COMPILER_NAME STREQUAL "nvcc"))
  set(BLA_VENDOR Intel10_64lp)
  find_package(BLAS)
  if(BLAS_FOUND) # in some systems with MKL, regular BLAS headers need to be found for it to work
    message("Multi/BLAS: MKL environment detected")
    add_definitions(-DRETURN_BY_STACK)
  else()
    message("Multi/BLAS: MKL environment not detected, looking for other BLAS")
    unset(BLA_VENDOR)
    find_package(BLAS REQUIRED)
    find_path(BLAS_INCLUDE_DIRS cblas.h
      /usr/include
      /usr/local/include
      $ENV{BLAS_HOME}/include
    )
  endif()

  find_path(BLAS_INCLUDE_DIRS cblas.h
    /usr/include
    /usr/local/include
    $ENV{BLAS_HOME}/include)

  link_libraries(${BLAS_LIBRARIES})
else()
  link_libraries("-lblas") # cmake cannot detect BLAS with pgi/nvc++ but it ships with its own version
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
  add_definitions(-DRETURN_BY_STACK)
endif()

include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(test)

