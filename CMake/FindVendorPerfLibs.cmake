# Copyright (c) 2026, The VendorPerfLibs Authors
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# FindVendorPerfLibs.cmake
#
# Searches for Vendor Performance Libraries (such as Intel MKL, AMD AOCL, or generic
# equivalents like Netlib LAPACK and FFTW) and provides a unified interface.
#
# This module supports the following components:
# - lapack : Linear Algebra PACKage
# - fft    : Fast Fourier Transform
# - vml    : Vector Math Library
#
# Example usage:
#   find_package(VendorPerfLibs COMPONENTS lapack fft vml REQUIRED)
#
# This module creates the following CMake imported targets (if their
# respective components are found):
# - VPL::lapack
# - VPL::fft
# - VPL::vml
#
# Input variables:
# - VPL_ID: Specifies the vendor performance library family to search for.
#   Acceptable values are:
#     - "IntelMKL" : Intel Math Kernel Library
#     - "AOCL"     : AMD Optimizing CPU Libraries
#     - "Generic"  : Generic libraries (e.g., standard BLAS/LAPACK and FFTW)
#   Note: If VPL_ID is not provided, the module will attempt to auto-detect
#   the appropriate vendor by checking for the presence of the MKL core or AOCL utils libraries.
#
# - VPL_OMP: If ON, OpenMP threading is requested. If OFF (default), sequential is used.
#   Note: If ON, you must call find_package(OpenMP) before finding VendorPerfLibs.
#   It is recommended to set the following in the project top-level for consistent selection of threading.
#   option(VPL_OMP "Use OpenMP threading for Vendor Performance Libraries" ${<project_OpenMP_variable>})
#
# Advanced Component-Specific Variables:
# - VPL_lapack_ID: Overrides VPL_ID specifically for the LAPACK component.
# - VPL_fft_ID: Overrides VPL_ID specifically for the FFT component.
# - VPL_vml_ID: Overrides VPL_ID specifically for the VML component.
#

include(FindPackageHandleStandardArgs)
include(CMakePushCheckState)

if(NOT (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED))
  message(FATAL_ERROR "VendorPerfLibs: C or CXX compiler must be loaded before calling find_package(VendorPerfLibs)")
endif()

cmake_push_check_state()

if(VPL_OMP)
  if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "VendorPerfLibs: VPL_OMP is ON but OpenMP_FOUND is false. Please invoke find_package(OpenMP) before VendorPerfLibs.")
  endif()
  if(CMAKE_Fortran_COMPILER_LOADED AND OpenMP_Fortran_FOUND)
    set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_Fortran_FLAGS})
  elseif(CMAKE_C_COMPILER_LOADED AND OpenMP_C_FOUND)
    set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_C_FLAGS})
  elseif(CMAKE_CXX_COMPILER_LOADED AND OpenMP_CXX_FOUND)
    set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_CXX_FLAGS})
  endif()
endif()

set(_VPL_VALID_IDS "IntelMKL" "AOCL" "Generic")
function(check_VPL_ID var_name id_to_check)
  if(NOT id_to_check IN_LIST _VPL_VALID_IDS)
    message(FATAL_ERROR "VendorPerfLibs: Unknown ${var_name} '${id_to_check}'. Acceptable values are: 'IntelMKL', 'Generic'")
  endif()
endfunction()

macro(speculateVendor)
  find_library(_MKL_CORE_TEST_LIB NAMES mkl_core
    HINTS
      "$ENV{MKLROOT}/lib/intel64"
  )

  find_library(AOCL_UTILS_LIB NAMES aoclutils)

  if(_MKL_CORE_TEST_LIB)
    set(VPL_ID_GUESS "IntelMKL")
  elseif(AOCL_UTILS_LIB)
    set(VPL_ID_GUESS "AOCL")
  else()
    set(VPL_ID_GUESS "Generic")
  endif()
endmacro()

if(NOT VPL_ID)
  speculateVendor()
else()
  check_VPL_ID("VPL_ID" "${VPL_ID}")
endif()

set(VPL_ID "${VPL_ID_GUESS}" CACHE STRING "Vendor Performance Library ID (IntelMKL, Generic)")

set(_find_package_args)
if(VendorPerfLibs_FIND_QUIETLY)
  list(APPEND _find_package_args QUIET)
endif()

function(warn_MKL_Fortran_ABI_issue CORE_LIB_VAR)
  # https://gitlab.kitware.com/cmake/cmake/-/merge_requests/12479
  # Switch to GNU Fortran ABI regarding how functions return complex numbers and how characters are passed (but not on Apple, where MKL does not provide it).
  # GNU and LLVMFlang families of compilers follow modern C99 _Complex ABI convention.
  # Intel, IntelLLVM, NVHPC follows the legacy convention.
  if(CMAKE_Fortran_COMPILER_LOADED AND CMAKE_VERSION VERSION_LESS "4.5")
    if(NOT (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") AND NOT APPLE)
      string(REPLACE "mkl_intel_lp64" "mkl_gf_lp64" CORE_LIBS "${${CORE_LIB_VAR}}")
    else()
      set(CORE_LIBS "${${CORE_LIB_VAR}}")
    endif()
    if(NOT VendorPerfLibs_FIND_QUIETLY AND NOT CORE_LIBS STREQUAL "${${CORE_LIB_VAR}}")
      message(WARNING "Potential incorrect selection of MKL ABI layer in variable ${CORE_LIB_VAR}.\n"
                      "    Original  : ${${CORE_LIB_VAR}}.\n"
                      "    Preferred : ${CORE_LIBS}\n"
                      "Use CMake >=4.5")
    endif()
  endif()
endfunction()

macro(find_VPL_core)
  set(VPL_CORE_FOUND TRUE)
  set(VPL_UTILS)
  if(VPL_ID STREQUAL "IntelMKL")
    # MKL core library support BLAS/LAPACK and FFT.
    # Thus we require VendorPerfLibs_INCLUDE_DIR and VPL_CORE_LIBRARIES being set
    list(APPEND _vpl_required_vars VendorPerfLibs_INCLUDE_DIR VPL_CORE_LIBRARIES)

    if(NOT VPL_OMP)
      set(BLA_VENDOR "Intel10_64lp_seq")
    else()
      set(BLA_VENDOR "Intel10_64lp")
    endif()
    find_package(BLAS ${_find_package_args})

    # Try to find MKL include directory
    find_path(VendorPerfLibs_INCLUDE_DIR NAMES mkl.h
      HINTS
        "$ENV{MKLROOT}/include"
      PATH_SUFFIXES mkl
    )

    if(BLAS_FOUND AND VendorPerfLibs_INCLUDE_DIR)
      warn_MKL_Fortran_ABI_issue(BLAS_LIBRARIES)
      set(VPL_CORE_LIBRARIES ${BLAS_LIBRARIES})
    else()
      set(VPL_CORE_FOUND FALSE)
      if(NOT VendorPerfLibs_FIND_QUIETLY)
        message(WARNING "Intel MKL not found. Please set the MKL root directory via CMake variable CMAKE_PREFIX_PATH or environment variable MKLROOT.")
      endif()
    endif()
  elseif(VPL_ID STREQUAL "AOCL")
    list(APPEND _vpl_required_vars VendorPerfLibs_INCLUDE_DIR)

    find_path(VendorPerfLibs_INCLUDE_DIR NAMES blis/blis.h blis.h)

    if(NOT VendorPerfLibs_INCLUDE_DIR)
      set(VPL_CORE_FOUND FALSE)
      if(NOT VendorPerfLibs_FIND_QUIETLY)
        message(WARNING "AMD AOCL include directory not found. Please set the AOCL root directory via CMake variable CMAKE_PREFIX_PATH.")
      endif()
    endif()

    find_library(AOCL_UTILS_LIB NAMES aoclutils)
    if(AOCL_UTILS_LIB)
      set(VPL_UTILS ${AOCL_UTILS_LIB})
    endif()
  endif()

  if(VPL_CORE_FOUND)
    list(APPEND _vpl_lib_found_ids "core(${VPL_ID})")
  endif()

  if(VPL_UTILS_LIB)
    list(APPEND _vpl_lib_found_ids "utils(${VPL_ID})")
  endif()
endmacro()

macro(find_VPL_lapack)
  set(VPL_lapack_ID ${VPL_ID} CACHE STRING "Vendor LAPACK ID (IntelMKL, Generic)")
  check_VPL_ID("VPL_lapack_ID" "${VPL_lapack_ID}")
  list(APPEND _vpl_required_vars LAPACK_LIBRARIES)

  if(VPL_lapack_ID STREQUAL "IntelMKL")
    if(NOT VPL_ID STREQUAL "IntelMKL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_lapack_ID is IntelMKL but VPL_ID is not IntelMKL. Unsupported.")
    endif()
    if(NOT VPL_OMP)
      set(BLA_VENDOR "Intel10_64lp_seq")
    else()
      set(BLA_VENDOR "Intel10_64lp")
    endif()
    find_package(LAPACK ${_find_package_args})
    if(LAPACK_FOUND)
      warn_MKL_Fortran_ABI_issue(LAPACK_LIBRARIES)
    endif()
  elseif(VPL_lapack_ID STREQUAL "AOCL")
    if(NOT VPL_ID STREQUAL "AOCL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_lapack_ID is AOCL but VPL_ID is not AOCL. Unsupported.")
    endif()
    if(NOT VPL_OMP)
      set(BLA_VENDOR "AOCL")
    else()
      set(BLA_VENDOR "AOCL_mt")
    endif()
    find_package(LAPACK ${_find_package_args})
  else()
    find_package(LAPACK ${_find_package_args})
  endif()

  if(LAPACK_FOUND)
    list(APPEND _vpl_lib_found_ids "lapack(${VPL_lapack_ID})")
    set(VendorPerfLibs_lapack_FOUND TRUE)
  else()
    set(VendorPerfLibs_lapack_FOUND FALSE)
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(WARNING "LAPACK for VPL_lapack_ID '${VPL_lapack_ID}' with OpenMP ${VPL_OMP}, not found")
    endif()
  endif()
endmacro()

macro(find_VPL_fft)
  set(VPL_fft_ID ${VPL_ID} CACHE STRING "Vendor FFT ID (IntelMKL, Generic)")
  check_VPL_ID("VPL_fft_ID" "${VPL_fft_ID}")
  list(APPEND _vpl_required_vars VendorPerfLibs_FFTW_INCLUDE_DIR VPL_FFT_LIBRARIES)

  set(VPL_FFT_FOUND TRUE)
  if(VPL_fft_ID STREQUAL "IntelMKL")
    if(NOT VPL_ID STREQUAL "IntelMKL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_fft_ID is IntelMKL but VPL_ID is not IntelMKL. Unsupported.")
    endif()
    find_path(VendorPerfLibs_FFTW_INCLUDE_DIR NAMES fftw3.f03
      HINTS
        "$ENV{MKLROOT}/include"
      PATHS
        /opt/intel/oneapi/mkl/latest/include
        /opt/intel/mkl/include
      PATH_SUFFIXES fftw mkl/fftw
    )
    set(VPL_FFT_LIBRARIES ${VPL_CORE_LIBRARIES})
    if(NOT VendorPerfLibs_FFTW_INCLUDE_DIR)
      set(VPL_FFT_FOUND FALSE)
    endif()
  else()
    # search for fftw
    if(VPL_OMP)
      find_package(FFTW ${_find_package_args} COMPONENTS seq omp)
      if(NOT FFTW_FOUND)
        find_package(FFTW ${_find_package_args} COMPONENTS seq)
      endif()
    else()
      find_package(FFTW ${_find_package_args} COMPONENTS seq)
    endif()
    if(FFTW_FOUND)
      set(VPL_FFT_FOUND TRUE)
      set(VendorPerfLibs_FFTW_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
      set(VPL_FFT_LIBRARIES ${FFTW_LIBRARIES})
    else()
      set(VPL_FFT_FOUND FALSE)
    endif()
  endif()

  if(VPL_FFT_FOUND)
    list(APPEND _vpl_lib_found_ids "fft(${VPL_fft_ID})")
    set(VendorPerfLibs_fft_FOUND TRUE)
  else()
    set(VendorPerfLibs_fft_FOUND FALSE)
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(WARNING "FFT for VPL_fft_ID '${VPL_fft_ID}' with OpenMP ${VPL_OMP}, not found")
    endif()
  endif()
endmacro()

macro(find_VPL_vml)
  set(VPL_vml_ID ${VPL_ID} CACHE STRING "Vendor VML ID (IntelMKL, Generic)")
  check_VPL_ID("VPL_vml_ID" "${VPL_vml_ID}")

  set(VPL_VML_FOUND TRUE)
  if(VPL_vml_ID STREQUAL "IntelMKL")
    if(NOT VPL_ID STREQUAL "IntelMKL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_vml_ID is IntelMKL but VPL_ID is not IntelMKL. Unsupported.")
    endif()
    # VML is part of MKL core. No additional libraries needed beyond core MKL libraries.
    set(VPL_VML_LIBRARIES ${VPL_CORE_LIBRARIES})
  elseif(VPL_vml_ID STREQUAL "AOCL")
    if(NOT VPL_ID STREQUAL "AOCL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_vml_ID is AOCL but VPL_ID is not AOCL. Unsupported.")
    endif()
    find_library(VPL_VML_LIBRARIES NAMES alm)
    if(NOT VPL_VML_LIBRARIES)
      set(VPL_VML_FOUND FALSE)
    endif()
  else()
    # Generic VML is not currently supported (e.g. no direct open source drop-in)
    set(VPL_VML_FOUND FALSE)
  endif()

  if(VPL_VML_FOUND)
    list(APPEND _vpl_lib_found_ids "vml(${VPL_vml_ID})")
    set(VendorPerfLibs_vml_FOUND TRUE)
  else()
    set(VendorPerfLibs_vml_FOUND FALSE)
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(WARNING "VML for VPL_vml_ID '${VPL_vml_ID}' not found or unsupported")
    endif()
  endif()
endmacro()

set(_vpl_lib_found_ids)
set(_vpl_required_vars _vpl_lib_found_ids)
find_VPL_core()

if(NOT VendorPerfLibs_FIND_COMPONENTS OR "lapack" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
  find_VPL_lapack()
endif()

if(NOT VendorPerfLibs_FIND_COMPONENTS OR "fft" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
  find_VPL_fft()
endif()

if(NOT VendorPerfLibs_FIND_COMPONENTS OR "vml" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
  find_VPL_vml()
endif()

find_package_handle_standard_args(VendorPerfLibs
  REQUIRED_VARS ${_vpl_required_vars}
  HANDLE_COMPONENTS
)

if(VendorPerfLibs_FOUND)
  # Create vpl_lapack
  if(NOT VendorPerfLibs_FIND_COMPONENTS OR "lapack" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
    if(NOT TARGET vpl_lapack)
      add_library(vpl_lapack INTERFACE IMPORTED)
      set_target_properties(vpl_lapack PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${VendorPerfLibs_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES};${VPL_UTILS}"
      )
      add_library(VPL::lapack ALIAS vpl_lapack)
    endif()
  endif()

  # Create vpl_fft
  if(NOT VendorPerfLibs_FIND_COMPONENTS OR "fft" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
    if(NOT TARGET vpl_fft)
      add_library(vpl_fft INTERFACE IMPORTED)
      set_target_properties(vpl_fft PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${VendorPerfLibs_INCLUDE_DIR};${VendorPerfLibs_FFTW_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${VPL_FFT_LIBRARIES};${VPL_UTILS}"
      )
      add_library(VPL::fft ALIAS vpl_fft)
    endif()
  endif()

  # Create vpl_vml
  if(NOT VendorPerfLibs_FIND_COMPONENTS OR "vml" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
    if(NOT TARGET vpl_vml)
      add_library(vpl_vml INTERFACE IMPORTED)
      set_target_properties(vpl_vml PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${VendorPerfLibs_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${VPL_VML_LIBRARIES};${VPL_UTILS}"
      )
      add_library(VPL::vml ALIAS vpl_vml)
    endif()
  endif()
endif()

cmake_pop_check_state()

# Scope cleanup
unset(_find_package_args)
unset(_vpl_required_vars)
unset(_vpl_lib_found_ids)
unset(_VPL_VALID_IDS)
