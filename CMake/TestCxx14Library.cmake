
# Test the compiler is configured with a C++14 library

try_compile(CXX14_LIBRARY_OKAY ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/CMake/try_cxx14_library.cpp)
message("C++ 14 library support = ${CXX14_LIBRARY_OKAY}")

if (NOT CXX14_LIBRARY_OKAY)
  message("C++ 14 library support not found")
  message("compiler is ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU") 
    message("Compiler detected is g++.  Use version 5.0 or newer for a C++14 compatible library")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("Compiler deteceted is clang++. If not using libcxx, ensure a g++ version greater than 5.0 is on the path")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("Compiler detected is icpc.  Ensure a g++ version greater than 5.0 is on the path.
                 Or use the --blah switch")
  endif()
  message(FATAL_ERROR "stopping")
endif()

