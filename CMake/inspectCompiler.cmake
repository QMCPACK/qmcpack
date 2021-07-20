# Map the compiler to the internal compiler flag/option customization
# COMPILER is defined upon return
#
# Note: vendor compilers can be just rebranded customized Clang compiler.
# It requires more recent CMake to handle it properly. We need to handle such cases for older CMake.

execute_process(
  COMMAND ${CMAKE_CXX_COMPILER} --version
  RESULT_VARIABLE VERSION_QUERY_RETURN
  OUTPUT_VARIABLE VERSION_QUERY_OUTPUT)

if(VERSION_QUERY_RETURN EQUAL 0)
  if(CMAKE_VERSION VERSION_LESS 3.20
     AND VERSION_QUERY_OUTPUT MATCHES "Intel"
     AND VERSION_QUERY_OUTPUT MATCHES "oneAPI")
    set(INTEL_ONEAPI_COMPILER_FOUND TRUE)
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(COMPILER GNU)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")
  set(COMPILER IBM)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR INTEL_ONEAPI_COMPILER_FOUND)
  set(COMPILER Intel)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "PGI" OR CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")
  set(COMPILER PGI)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
  set(COMPILER Cray)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(COMPILER Clang)
else()
  message("${CMAKE_CXX_COMPILER_ID}")
  message(WARNING "Unknown C/C++ compiler ${CMAKE_CXX_COMPILER_ID}, default flags will be used")
endif()
message(STATUS "C++ Compiler is identified by QMCPACK as : ${COMPILER}")

#------------------------------------
# Include compiler-specific cmake file
#------------------------------------
if(COMPILER MATCHES "IBM")
  include(${PROJECT_CMAKE}/IBMCompilers.cmake)
elseif(COMPILER MATCHES "Intel")
  include(${PROJECT_CMAKE}/IntelCompilers.cmake)
elseif(COMPILER MATCHES "GNU")
  include(${PROJECT_CMAKE}/GNUCompilers.cmake)
elseif(COMPILER MATCHES "Clang")
  include(${PROJECT_CMAKE}/ClangCompilers.cmake)
elseif(COMPILER MATCHES "PGI")
  include(${PROJECT_CMAKE}/PGICompilers.cmake)
else()
  message(WARNING "No default file for compiler (${COMPILER})")
endif()
