# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
IntelDPCPPConfig
-------

DPCPP Library to verify DPCPP/SYCL compatability of CMAKE_CXX_COMPILER
and passes relevant compiler flags.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``IntelDPCPP_FOUND``
  True if the system has the DPCPP library.
``SYCL_LANGUAGE_VERSION``
  The SYCL language spec version by Compiler.
``SYCL_INCLUDE_DIR``
  Include directories needed to use SYCL.
``SYCL_IMPLEMENTATION_ID``
  The SYCL compiler variant.
``SYCL_FLAGS``
  SYCL specific flags for the compiler.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``SYCL_INCLUDE_DIR``
  The directory containing ``sycl.hpp``.
``SYCL_LIBRARY_DIR``
  The path to the SYCL library.
``SYCL_FLAGS``
  SYCL specific flags for the compiler.
``SYCL_LANGUAGE_VERSION``
  The SYCL language spec version by Compiler.


.. note::

  For now, user needs to set -DCMAKE_CXX_COMPILER or environment of
  CXX pointing to SYCL compatible compiler  ( eg: icx, clang++, icpx)

  Note: do not set to DPCPP compiler. If set to a Compiler family
  that supports dpcpp ( eg: IntelLLVM) both DPCPP and SYCL
  features are enabled.

  And add this package to user's Cmake config file.

  .. code-block:: cmake

    find_package(IntelDPCPP REQUIRED)

#]=======================================================================]

include(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
  # TODO add dependency package module checks, if any
endif()


# TODO: can't use find_program to override the CMAKE_CXX_COMPILER as
# Platform/ files are executed, potentially for a different compiler.
# Safer approach is to make user to define CMAKE_CXX_COMPILER.

string(COMPARE EQUAL "${CMAKE_CXX_COMPILER}" "" nocmplr)
if(nocmplr)
  set(IntelDPCPP_FOUND False)
  set(SYCL_REASON_FAILURE "SYCL: CMAKE_CXX_COMPILER not set!!")
  set(IntelDPCPP_NOT_FOUND_MESSAGE "${SYCL_REASON_FAILURE}")
endif()

# Check for known compiler family that supports SYCL

if( NOT "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xClang" AND
    NOT "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntelLLVM")
   set(IntelDPCPP_FOUND False)
   set(SYCL_REASON_FAILURE "Unsupported compiler family ${CMAKE_CXX_COMPILER_ID} and compiler ${CMAKE_CXX_COMPILER}!!")
   set(IntelDPCPP_NOT_FOUND_MESSAGE "${SYCL_REASON_FAILURE}")
   return()
endif()

# Assume that CXX Compiler supports SYCL and then test to verify.
set(SYCL_COMPILER ${CMAKE_CXX_COMPILER})


# Function to write a test case to verify SYCL features.

function(SYCL_FEATURE_TEST_WRITE src)

  set(pp_if "#if")
  set(pp_endif "#endif")

  set(SYCL_TEST_CONTENT "")
  string(APPEND SYCL_TEST_CONTENT "#include <iostream>\nusing namespace std;\n")
  string(APPEND SYCL_TEST_CONTENT "int main(){\n")

  # Feature tests goes here

  string(APPEND SYCL_TEST_CONTENT "${pp_if} defined(SYCL_LANGUAGE_VERSION)\n")
  string(APPEND SYCL_TEST_CONTENT "cout << \"SYCL_LANGUAGE_VERSION=\"<<SYCL_LANGUAGE_VERSION<<endl;\n")
  string(APPEND SYCL_TEST_CONTENT "${pp_endif}\n")

  string(APPEND SYCL_TEST_CONTENT "return 0;}\n")

  file(WRITE ${src} "${SYCL_TEST_CONTENT}")

endfunction()

# Function to Build the feature check test case.

function(SYCL_FEATURE_TEST_BUILD TEST_SRC_FILE TEST_EXE)

  # Convert CXX Flag string to list
  set(SYCL_CXX_FLAGS_LIST "${SYCL_CXX_FLAGS}")
  separate_arguments(SYCL_CXX_FLAGS_LIST)

  # Spawn a process to build the test case.
  execute_process(
    COMMAND "${SYCL_COMPILER}"
    ${SYCL_CXX_FLAGS_LIST}
    ${TEST_SRC_FILE}
    "-o"
    ${TEST_EXE}
    WORKING_DIRECTORY ${SYCL_TEST_DIR}
    OUTPUT_VARIABLE output ERROR_VARIABLE output
    OUTPUT_FILE ${SYCL_TEST_DIR}/Compile.log
    RESULT_VARIABLE result
    TIMEOUT 20
    )

  # Verify if test case build properly.
  if(result)
    message("SYCL feature test compile failed!")
    message("compile output is: ${output}")
  endif()

  # TODO: what to do if it doesn't build

endfunction()

# Function to run the test case to generate feature info.

function(SYCL_FEATURE_TEST_RUN TEST_EXE)

  # Spawn a process to run the test case.

  execute_process(
    COMMAND ${TEST_EXE}
    WORKING_DIRECTORY ${SYCL_TEST_DIR}
    OUTPUT_VARIABLE output ERROR_VARIABLE output
    RESULT_VARIABLE result
    TIMEOUT 20
    )

  # Verify the test execution output.
  if(test_result)
    set(IntelDPCPP_FOUND False)
    set(SYCL_REASON_FAILURE "SYCL: feature test execution failed!!")
  endif()
  # TODO: what iff the result is false.. error or ignore?

  set( test_result "${result}" PARENT_SCOPE)
  set( test_output "${output}" PARENT_SCOPE)

endfunction()


# Function to extract the information from test execution.
function(SYCL_FEATURE_TEST_EXTRACT test_output)

  string(REGEX REPLACE "\n" ";" test_output_list "${test_output}")

  set(SYCL_LANGUAGE_VERSION "")
  foreach(strl ${test_output_list})
     if(${strl} MATCHES "^SYCL_LANGUAGE_VERSION=([A-Za-z0-9_]+)$")
       string(REGEX REPLACE "^SYCL_LANGUAGE_VERSION=" "" extracted_sycl_lang "${strl}")
       set(SYCL_LANGUAGE_VERSION ${extracted_sycl_lang})
     endif()
  endforeach()

  set(SYCL_LANGUAGE_VERSION "${SYCL_LANGUAGE_VERSION}" PARENT_SCOPE)
endfunction()

if(SYCL_COMPILER)
  # TODO ensure CMAKE_LINKER and CMAKE_CXX_COMPILER are same/supports SYCL.
  # set(CMAKE_LINKER ${SYCL_COMPILER})

  # use REALPATH to resolve symlinks
  get_filename_component(_REALPATH_SYCL_COMPILER "${SYCL_COMPILER}" REALPATH)
  get_filename_component(SYCL_BIN_DIR "${_REALPATH_SYCL_COMPILER}" DIRECTORY)
  get_filename_component(SYCL_PACKAGE_DIR "${SYCL_BIN_DIR}" DIRECTORY CACHE)

  # Find Include path from binary
  find_file(SYCL_INCLUDE_DIR
    NAMES
      include
    HINTS
      ${SYCL_PACKAGE_DIR}
    NO_DEFAULT_PATH
      )

  # Find Library directory
  find_file(SYCL_LIBRARY_DIR
    NAMES
      lib lib64
    HINTS
      ${SYCL_PACKAGE_DIR}
    NO_DEFAULT_PATH
      )

endif()


set(SYCL_FLAGS "")

# Based on Compiler ID, add support for SYCL
if( "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xClang" OR
    "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntelLLVM")
  set(SYCL_FLAGS "-fsycl")
endif()

# Based on Compiler ID, add support for DPCPP
if( "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntelLLVM")
  list(PREPEND SYCL_FLAGS "--dpcpp")
endif()

# TODO verify if this is needed
# Windows: Add Exception handling
if(WIN32)
  list(APPEND SYCL_FLAGS "/EHsc")
endif()

set(SYCL_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYCL_FLAGS}")

# And now test the assumptions.

# Create a clean working directory.
set(SYCL_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/TESTSYCL")
file(REMOVE_RECURSE ${SYCL_TEST_DIR})
file(MAKE_DIRECTORY ${SYCL_TEST_DIR})

# Create the test source file
set(TEST_SRC_FILE "${SYCL_TEST_DIR}/sycl_features.cpp")
set(TEST_EXE "${TEST_SRC_FILE}.exe")
SYCL_FEATURE_TEST_WRITE(${TEST_SRC_FILE})

# Build the test and create test executable
SYCL_FEATURE_TEST_BUILD(${TEST_SRC_FILE} ${TEST_EXE})

# Execute the test to extract information
SYCL_FEATURE_TEST_RUN(${TEST_EXE})

# Extract test output for information
SYCL_FEATURE_TEST_EXTRACT(${test_output})

# As per specification, all the SYCL compatible compilers should
# define macro  SYCL_LANGUAGE_VERSION
string(COMPARE EQUAL "${SYCL_LANGUAGE_VERSION}" "" nosycllang)
if(nosycllang)
  set(IntelDPCPP_FOUND False)
  set(SYCL_REASON_FAILURE "SYCL: It appears that the ${CMAKE_CXX_COMPILER} does not support SYCL")
  set(IntelDPCPP_NOT_FOUND_MESSAGE "${SYCL_REASON_FAILURE}")
endif()

# Placeholder for identifying various implemenations of SYCL compilers.
# for now, set to the CMAKE_CXX_COMPILER_ID
set(SYCL_IMPLEMENTATION_ID "${CMAKE_CXX_COMPILER_ID}")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYCL_FLAGS}")
find_library(
  SYCL_LIBRARY
  NAMES
    sycl
  HINTS
    ${SYCL_LIBRARY_DIR}
  NO_DEFAULT_PATH
)


message(STATUS "The SYCL compiler is ${SYCL_COMPILER}")
message(STATUS "The SYCL Flags are ${SYCL_FLAGS}")
message(STATUS "The SYCL Language Version is ${SYCL_LANGUAGE_VERSION}")

find_package_handle_standard_args(
  IntelDPCPP
  FOUND_VAR IntelDPCPP_FOUND
  REQUIRED_VARS SYCL_INCLUDE_DIR SYCL_LIBRARY_DIR SYCL_FLAGS SYCL_LIBRARY
  VERSION_VAR SYCL_LANGUAGE_VERSION
  REASON_FAILURE_MESSAGE "${SYCL_REASON_FAILURE}")

# Include in Cache
set(SYCL_LANGUAGE_VERSION "${SYCL_LANGUAGE_VERSION}" CACHE STRING "SYCL Language version")
set(SYCL_INCLUDE_DIR "${SYCL_INCLUDE_DIR}" CACHE FILEPATH "SYCL Include directory")
set(SYCL_LIBRARY_DIR "${SYCL_LIBRARY_DIR}" CACHE FILEPATH "SYCL Library Directory")
set(SYCL_FLAGS "${SYCL_FLAGS}" CACHE STRING "SYCL flags for the compiler")

# sycl target for device compilation
add_library(OneAPI::DPCPP-device INTERFACE IMPORTED)
target_compile_options(OneAPI::DPCPP-device INTERFACE ${SYCL_FLAGS})
target_link_options(OneAPI::DPCPP-device INTERFACE ${SYCL_FLAGS})

# sycl target for host APIs
add_library(OneAPI::DPCPP-host INTERFACE IMPORTED)
target_include_directories(OneAPI::DPCPP-host INTERFACE ${SYCL_INCLUDE_DIR})
target_link_libraries(OneAPI::DPCPP-host INTERFACE ${SYCL_LIBRARY})
