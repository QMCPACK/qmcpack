# Use cmake/ctest for building and testing boost mpi3 tests
cmake_minimum_required(VERSION 3.18.4) # 3.10 for FindMPI, 3.12 for MPI_ROOT
project(boostmpi3nccl_tests VERSION 0.1 LANGUAGES CXX)

enable_testing()
include(CTest)

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
link_libraries(${MPI_CXX_LIBRARIES})

if (NOT MPIEXEC)
  set(MPIEXEC mpirun)
endif()


# Set BOOST_ROOT if boost is in a non-standard location
find_package(Boost REQUIRED)
if(Boost_FOUND)
  set(HAVE_LIBBOOST 1)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  message(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
endif()

#if(ENABLE_CUDA)
	enable_language(CUDA)
	find_package(CUDA 11.7.0 REQUIRED)
#	set(CXX_COMPILER ${CMAKE_CUDA_COMPILER})
#	set(CXX_FLAGS ${CMAKE_CUDA_FLAGS})
#	string(APPEND CMAKE_CUDA_FLAGS " --forward-unknown-to-host-linker") # ${ARCH_FLAGS} -std=c++17 --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"")
#	add_compile_options(--forward-unknown-to-host-linker)
#endif()

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--display_error_number\" --std=c++17")

set(NCCL_INCLUDE_DIR $ENV{NCCL_INCLUDE_DIR} CACHE PATH "Folder contains NVIDIA NCCL headers")
set(NCCL_LIB_DIR $ENV{NCCL_LIB_DIR} CACHE PATH "Folder contains NVIDIA NCCL libraries")
set(NCCL_VERSION $ENV{NCCL_VERSION} CACHE STRING "Version of NCCL to build with")

if ($ENV{NCCL_ROOT_DIR})
  message(WARNING "NCCL_ROOT_DIR is deprecated. Please set NCCL_ROOT instead.")
endif()
list(APPEND NCCL_ROOT $ENV{NCCL_ROOT_DIR} ${CUDA_TOOLKIT_ROOT_DIR})
# Compatible layer for CMake <3.12. NCCL_ROOT will be accounted in for searching paths and libraries for CMake >=3.12.
list(APPEND CMAKE_PREFIX_PATH ${NCCL_ROOT})

find_path(NCCL_INCLUDE_DIRS
  NAMES nccl.h
  HINTS ${NCCL_INCLUDE_DIR})

if (USE_STATIC_NCCL)
  MESSAGE(STATUS "USE_STATIC_NCCL is set. Linking with static NCCL library.")
  SET(NCCL_LIBNAME "nccl_static")
  if (NCCL_VERSION)  # Prefer the versioned library if a specific NCCL version is specified
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a.${NCCL_VERSION}" ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
else()
  SET(NCCL_LIBNAME "nccl")
  if (NCCL_VERSION)  # Prefer the versioned library if a specific NCCL version is specified
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".so.${NCCL_VERSION}" ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
endif()

find_library(NCCL_LIBRARIES
  NAMES ${NCCL_LIBNAME}
  HINTS ${NCCL_LIB_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NCCL DEFAULT_MSG NCCL_INCLUDE_DIRS NCCL_LIBRARIES)

if(NCCL_FOUND)  # obtaining NCCL version and some sanity checks
  set (NCCL_HEADER_FILE "${NCCL_INCLUDE_DIRS}/nccl.h")
  message (STATUS "Determining NCCL version from ${NCCL_HEADER_FILE}...")
  set (OLD_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
  list (APPEND CMAKE_REQUIRED_INCLUDES ${NCCL_INCLUDE_DIRS})
  include(CheckCXXSymbolExists)
  check_cxx_symbol_exists(NCCL_VERSION_CODE nccl.h NCCL_VERSION_DEFINED)

  if (NCCL_VERSION_DEFINED)
    set(file "${PROJECT_BINARY_DIR}/detect_nccl_version.cc")
    file(WRITE ${file} "
      #include <iostream>
      #include <nccl.h>
      int main()
      {
        std::cout << NCCL_MAJOR << '.' << NCCL_MINOR << '.' << NCCL_PATCH << std::endl;
        int x;
        ncclGetVersion(&x);
        return x == NCCL_VERSION_CODE;
      }
")
    try_run(NCCL_VERSION_MATCHED compile_result ${PROJECT_BINARY_DIR} ${file}
          RUN_OUTPUT_VARIABLE NCCL_VERSION_FROM_HEADER
          CMAKE_FLAGS  "-DINCLUDE_DIRECTORIES=${NCCL_INCLUDE_DIRS}"
          LINK_LIBRARIES ${NCCL_LIBRARIES})
    if (NOT NCCL_VERSION_MATCHED)
      message(FATAL_ERROR "Found NCCL header version and library version do not match! \
(include: ${NCCL_INCLUDE_DIRS}, library: ${NCCL_LIBRARIES}) Please set NCCL_INCLUDE_DIR and NCCL_LIB_DIR manually.")
    endif()
    message(STATUS "NCCL version: ${NCCL_VERSION_FROM_HEADER}")
  else()
    message(STATUS "NCCL version < 2.3.5-5")
  endif ()
  set (CMAKE_REQUIRED_INCLUDES ${OLD_CMAKE_REQUIRED_INCLUDES})

  message(STATUS "Found NCCL (include: ${NCCL_INCLUDE_DIRS}, library: ${NCCL_LIBRARIES})")
  mark_as_advanced(NCCL_ROOT_DIR NCCL_INCLUDE_DIRS NCCL_LIBRARIES)
endif()

# This list is only a subset of files.  It might be useful to divide some of the tests into
#  different categories (MPI-1, MPI-2, etc., or something else)
set(TEST_SRCS
	nccl_constructor.cu
)

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}.x")
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_compile_features(${TEST_EXE} PRIVATE cxx_std_17)

#  if(ENABLE_CUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
#  endif()

  target_compile_options(${TEST_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Werror>) # -Wpedantic -Wunused -Wnon-virtual-dtor -Woverloaded-virtual -Wcast-qual -Wformat=2>) #-Wsign-conversion -Wshadow -Wconversion
  target_compile_options(${TEST_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Wall -Wextra>)
  target_include_directories(${TEST_EXE} PUBLIC "../../../..")

  target_include_directories(${TEST_EXE} SYSTEM PRIVATE /home/correaa/prj/alf/boost/multi/include)

  list(FIND NEED_BOOST_SERIALIZATION_SRCS ${TEST_FILE} NEED_BOOST_SERIALIZATION)
  if (NOT (${NEED_BOOST_SERIALIZATION} EQUAL -1))
    target_compile_options(${TEST_EXE} PUBLIC -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY)
#    target_link_libraries(${TEST_EXE} boost_serialization)
  endif()

  target_link_libraries(${TEST_EXE} pthread)

  target_include_directories(${TEST_EXE} PUBLIC ${NCCL_INCLUDE_DIRS})
  target_link_libraries(${TEST_EXE} ${NCCL_LIBRARIES})

  add_test(NAME ${TEST_EXE} COMMAND ${MPIEXEC} $ENV{MPI_OVERSUBSCRIBE} -n ${NPROC} ./${TEST_EXE})
  set_tests_properties(${TEST_EXE} PROPERTIES TIMEOUT 30)
endforeach()
