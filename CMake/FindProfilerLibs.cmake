# CMake note - complex conditionals in cmake_dependent_option must have spaces around parentheses
cmake_dependent_option(USE_NVTX_API "Enable/disable NVTX regions in CUDA code." OFF "ENABLE_TIMERS AND ENABLE_CUDA" OFF)
option(USE_VTUNE_API "Enable use of VTune ittnotify APIs" OFF)
cmake_dependent_option(USE_VTUNE_TASKS "USE VTune ittnotify task annotation" OFF "ENABLE_TIMERS AND USE_VTUNE_API" OFF)
option(USE_HPCTOOLKIT_API "Enable use of HPCToolkit start/stop APIs" OFF)

#-------------------------------------------------------------------
# set up NVTX library
#-------------------------------------------------------------------
if(USE_NVTX_API AND NOT QMC_CUDA2HIP)
  message(STATUS "Enabling use of CUDA NVTX APIs")

  set(_NVTX_FOUND FALSE)

  # 1) trying to find NVTX v3 (header-only)
  if(TARGET CUDA::nvtx3)
    message(STATUS "Using CUDA::nvtx3 (header-only NVTX v3)")
    link_libraries(CUDA::nvtx3)
    set(_NVTX_FOUND TRUE)
  endif()

  # 2) if CUDA::nvtx3 fails, try to find the nvtx3 headers manually.
  if(NOT _NVTX_FOUND)
    find_path(NVTX3_INCLUDE_DIR
      NAMES
        nvtx3/nvtx3.hpp
        nvtx3/nvToolsExt.h
      HINTS
        ${CUDA_TOOLKIT_ROOT_DIR}/include
        $ENV{CUDA_PATH}/include
    )
    if(NVTX3_INCLUDE_DIR)
      message(STATUS "Found NVTX v3 headers at: ${NVTX3_INCLUDE_DIR} (header-only)")
      include_directories(${NVTX3_INCLUDE_DIR})
      set(_NVTX_FOUND TRUE)
    endif()
  endif()

  # generate shim header for NVTX v3, if found
  if(TARGET CUDA::nvtx3 OR NVTX3_INCLUDE_DIR)
    set(NVTX_SHIM_DIR "${CMAKE_BINARY_DIR}/nvtx_compat")
    file(MAKE_DIRECTORY "${NVTX_SHIM_DIR}")
    file(WRITE "${NVTX_SHIM_DIR}/nvToolsExt.h"
"#pragma once
#include <nvtx3/nvToolsExt.h>"
    )
    include_directories(BEFORE "${NVTX_SHIM_DIR}")
    message(STATUS "NVTX: using v3 compat shim at ${NVTX_SHIM_DIR}/nvToolsExt.h")
  endif()

  # 3) fallback to legacy NVTX v2
  if(NOT _NVTX_FOUND)
    find_library(NVTX2_LIB
      NAME nvToolsExt
      HINTS
        ${CUDA_TOOLKIT_ROOT_DIR}
        $ENV{CUDA_PATH}
      PATH_SUFFIXES
        lib
        lib64
    )
    if(NVTX2_LIB)
      message(STATUS "Using legacy NVTX v2: ${NVTX2_LIB}")
      link_libraries(${NVTX2_LIB})
      set(_NVTX_FOUND TRUE)
    endif()
  endif()

  if(NOT _NVTX_FOUND)
    message(FATAL_ERROR "USE_NVTX_API set but NVTX not found")
  endif(NOT _NVTX_FOUND)
else()
  message(STATUS "CUDA NVTX APIs disabled")
endif(USE_NVTX_API AND NOT QMC_CUDA2HIP)

#-------------------------------------------------------------------
# set up ROCTX library
#-------------------------------------------------------------------
if(USE_NVTX_API AND QMC_CUDA2HIP)
  message(STATUS "Enabling use of ROCm rocTX APIs")
  find_library(
    ROCTX_API_LIB NAME libroctx64.so
    PATH_SUFFIXES lib lib64)
  if(NOT ROCTX_API_LIB)
    message(FATAL_ERROR "USE_ROCTX_API set but ROCTX_API_LIB not found")
  endif(NOT ROCTX_API_LIB)
  message("ROCm rocTX library: ${ROCTX_API_LIB}")
  link_libraries(${ROCTX_API_LIB})
  else()
    message(STATUS "ROCm rocTX APIs disabled")
endif(USE_NVTX_API AND QMC_CUDA2HIP)

#-------------------------------------------------------------------
# set up VTune ittnotify library
#-------------------------------------------------------------------
if(USE_VTUNE_API)
  message(STATUS "Enabling use of VTune ittnotify APIs")
  find_path(
    VTUNE_ITTNOTIFY_INCLUDE_DIR ittnotify.h
    HINTS ${VTUNE_ROOT} $ENV{VTUNE_ROOT}
    PATH_SUFFIXES include REQUIRED)
  message(STATUS "Found VTUNE_ITTNOTIFY_INCLUDE_DIR ${VTUNE_ITTNOTIFY_INCLUDE_DIR}")
  find_library(
    VTUNE_ITTNOTIFY_LIBRARY ittnotify
    HINTS ${VTUNE_ROOT} $ENV{VTUNE_ROOT}
    PATH_SUFFIXES lib64 lib REQUIRED)
  message(STATUS "Found VTUNE_ITTNOTIFY_LIBRARY ${VTUNE_ITTNOTIFY_LIBRARY}")

  include_directories(${VTUNE_ITTNOTIFY_INCLUDE_DIR})
  link_libraries(${VTUNE_ITTNOTIFY_LIBRARY})
  if(USE_VTUNE_TASKS)
    message(STATUS "VTune ittnotify tasks enabled")
  endif()
else()
  message(STATUS "VTune ittnotify APIs disabled")
endif()

#-------------------------------------------------------------------
# set up HPCToolkit library
#-------------------------------------------------------------------
if(USE_HPCTOOLKIT_API)
  message(STATUS "Enabling use of HPCToolkit start/stop APIs")
  find_path(
    HPCTOOLKIT_INCLUDE_DIR hpctoolkit.h
    HINTS ${HPCTOOLKIT_ROOT} $ENV{HPCTOOLKIT_ROOT}
    PATH_SUFFIXES include REQUIRED)
  message(STATUS "Found HPCTOOLKIT_INCLUDE_DIR ${HPCTOOLKIT_INCLUDE_DIR}")
  find_library(
    HPCTOOLKIT_LIBRARY hpctoolkit
    HINTS ${HPCTOOLKIT_ROOT} $ENV{HPCTOOLKIT_ROOT}
    PATH_SUFFIXES lib64 lib REQUIRED)
  message(STATUS "Found HPCTOOLKIT_LIBRARY ${HPCTOOLKIT_LIBRARY}")

  include_directories(${HPCTOOLKIT_INCLUDE_DIR})
  link_libraries(${HPCTOOLKIT_LIBRARY})
else()
  message(STATUS "HPCToolkit start/stop APIs disabled")
endif()
