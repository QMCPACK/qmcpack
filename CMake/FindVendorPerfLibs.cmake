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
# This module is maintained at https://github.com/ye-luo/cmake_repo/tree/master/VendorPerfLibs
#
# Searches for Vendor Performance Libraries (such as Intel MKL, AMD AOCL, or generic
# equivalents like Netlib LAPACK and FFTW) and provides a unified interface.
#
# This module supports the following components:
# - lapack : Linear Algebra PACKage
# - fft    : Fast Fourier Transform
# - vml    : Vector Math Library
#
# Requirements:
# - The fft component for the Generic and AOCL vendors relies on a project-supplied FindFFTW.cmake
#   (it is not part of CMake). It must be reachable via CMAKE_MODULE_PATH and support the COMPONENTS
#   "seq" and "omp", and define FFTW_FOUND, FFTW_INCLUDE_DIR and FFTW_LIBRARIES.
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
# These imported targets are the supported interface of this module. Link against them rather than
# relying on any other variable set by this module.
#
# Result variables:
# - VendorPerfLibs_FOUND         : TRUE if the requested components were found.
# - VendorPerfLibs_<comp>_FOUND  : TRUE if component <comp> (lapack, fft or vml) was found.
#
# Input variables:
# - VPL_ID: Specifies the vendor performance library family to search for.
#   Acceptable values are:
#     - "IntelMKL"   : Intel Math Kernel Library
#     - "AOCL"       : AMD Optimizing CPU Libraries
#     - "ARMPL"      : Arm Performance Libraries
#     - "NVPL"       : NVIDIA Performance Libraries
#     - "CrayLibSci" : Cray Scientific Library
#     - "Generic"    : Generic libraries (e.g., standard BLAS/LAPACK and FFTW)
#   Note: If VPL_ID is not set, the following rules determine its value
#   (once a rule applies, the subsequent rules are ignored):
#     1. "Generic" if BLA_VENDOR was specified.
#     2. "IntelMKL" if the Intel MKL mkl_core library file was found.
#     3. "AOCL" if the AOCL aoclutils library file was found.
#     4. "ARMPL" if the ArmPL library file was found.
#     5. "NVPL" if the nvpl CMake package was found.
#     6. "CrayLibSci" if the Cray libsci pkg-config package was found.
#     7. "Generic" if all previous rules fail.
#
#   Note: Autodetection skips (with a warning) AOCL if CMake < 3.27 and NVPL if CMake < 4.1, since the
#   corresponding BLA_VENDOR values are unavailable. Setting VPL_ID=AOCL or VPL_ID=NVPL explicitly with an
#   older CMake is a fatal error.
#
# - VPL_OMP: If ON, OpenMP threading is requested. If OFF (default), sequential is used.
#   Note: If ON, you must call find_package(OpenMP) before finding VendorPerfLibs.
#   It is recommended to set the following in the project top-level for consistent selection of threading.
#   option(VPL_OMP "Use OpenMP threading for Vendor Performance Libraries" ${<project_OpenMP_variable>})
#
# Advanced Component-Specific Variables:
# - VPL_fft_ID: Overrides VPL_ID specifically for the FFT component.
# - VPL_vml_ID: Overrides VPL_ID specifically for the VML component.
#

include(FindPackageHandleStandardArgs)
include(FindPackageMessage)
include(CMakePushCheckState)
find_package(PkgConfig QUIET)

# Cray prefixes must be on CMAKE_PREFIX_PATH for FindBLAS/FindLAPACK (pkg-config) and pkg_check_modules.
# Save the caller's value, so that it can be restored if the libraries are not found (see end of this module).
set(_vpl_saved_CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}")
if(DEFINED ENV{CRAY_LIBSCI_PREFIX_DIR} AND NOT "$ENV{CRAY_LIBSCI_PREFIX_DIR}" IN_LIST CMAKE_PREFIX_PATH)
  list(APPEND CMAKE_PREFIX_PATH "$ENV{CRAY_LIBSCI_PREFIX_DIR}")
endif()

if(NOT (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED))
  message(FATAL_ERROR "VendorPerfLibs: C or CXX compiler must be loaded before calling find_package(VendorPerfLibs)")
endif()

cmake_push_check_state()

if(VPL_OMP)
  if(NOT OpenMP_FOUND)
    message(FATAL_ERROR "VendorPerfLibs: VPL_OMP is ON but OpenMP_FOUND is false. Please invoke find_package(OpenMP) before VendorPerfLibs.")
  endif()
  # CMAKE_REQUIRED_LINK_OPTIONS expects a list, whereas OpenMP_<lang>_FLAGS is a space-separated string.
  # (CMAKE_REQUIRED_LIBRARIES cannot be used here since FindBLAS/FindLAPACK overwrite it.)
  set(_vpl_omp_flags)
  if(CMAKE_Fortran_COMPILER_LOADED AND OpenMP_Fortran_FOUND)
    separate_arguments(_vpl_omp_flags NATIVE_COMMAND "${OpenMP_Fortran_FLAGS}")
  elseif(CMAKE_C_COMPILER_LOADED AND OpenMP_C_FOUND)
    separate_arguments(_vpl_omp_flags NATIVE_COMMAND "${OpenMP_C_FLAGS}")
  elseif(CMAKE_CXX_COMPILER_LOADED AND OpenMP_CXX_FOUND)
    separate_arguments(_vpl_omp_flags NATIVE_COMMAND "${OpenMP_CXX_FLAGS}")
  endif()
  set(CMAKE_REQUIRED_LINK_OPTIONS ${_vpl_omp_flags})
  unset(_vpl_omp_flags)
endif()

set(_VPL_VALID_IDS "IntelMKL" "AOCL" "ARMPL" "NVPL" "CrayLibSci" "Generic")
function(_vpl_check_ID var_name id_to_check)
  if(NOT id_to_check IN_LIST _VPL_VALID_IDS)
    message(FATAL_ERROR "VendorPerfLibs: Unknown ${var_name} '${id_to_check}'. Acceptable values are: 'IntelMKL', 'AOCL', 'ARMPL', 'NVPL', 'CrayLibSci', 'Generic'")
  endif()
endfunction()

# FindBLAS/FindLAPACK are driven through BLA_* variables, which this module sets.
# Save the caller's values on entry so that they can be restored for a Generic search, or on module exit
# when the vendor libraries were not found. When found, the BLA_* variables intentionally persist.
set(_VPL_BLA_VARS BLA_VENDOR BLA_THREAD BLA_PREFER_PKGCONFIG BLA_PKGCONFIG_BLAS BLA_PKGCONFIG_LAPACK)
macro(_vpl_save_BLA_vars)
  foreach(_vpl_bla_var IN LISTS _VPL_BLA_VARS)
    if(DEFINED ${_vpl_bla_var})
      set(_vpl_saved_${_vpl_bla_var} "${${_vpl_bla_var}}")
      set(_vpl_had_${_vpl_bla_var} TRUE)
    else()
      set(_vpl_had_${_vpl_bla_var} FALSE)
    endif()
  endforeach()
  unset(_vpl_bla_var)
endmacro()
macro(_vpl_restore_BLA_vars)
  foreach(_vpl_bla_var IN LISTS _VPL_BLA_VARS)
    if(_vpl_had_${_vpl_bla_var})
      set(${_vpl_bla_var} "${_vpl_saved_${_vpl_bla_var}}")
    else()
      unset(${_vpl_bla_var})
    endif()
  endforeach()
  unset(_vpl_bla_var)
endmacro()
macro(_vpl_clear_BLA_vars)
  foreach(_vpl_bla_var IN LISTS _VPL_BLA_VARS)
    unset(_vpl_saved_${_vpl_bla_var})
    unset(_vpl_had_${_vpl_bla_var})
  endforeach()
  unset(_vpl_bla_var)
endmacro()

_vpl_save_BLA_vars()

# MKL search hints from MKLROOT. Only add them when MKLROOT is set: otherwise "$ENV{MKLROOT}/include" would
# expand to "/include" and search the filesystem root.
set(_VPL_MKL_LIB_HINTS)
set(_VPL_MKL_INC_HINTS)
if(DEFINED ENV{MKLROOT} AND NOT "$ENV{MKLROOT}" STREQUAL "")
  set(_VPL_MKL_LIB_HINTS "$ENV{MKLROOT}/lib/intel64" "$ENV{MKLROOT}/lib")
  set(_VPL_MKL_INC_HINTS "$ENV{MKLROOT}/include")
endif()

macro(_vpl_speculate_vendor)
  find_library(_MKL_CORE_TEST_LIB NAMES mkl_core
    HINTS
      ${_VPL_MKL_LIB_HINTS}
  )

  find_library(_AOCL_UTILS_LIB NAMES aoclutils)

  find_library(_ARMPL_TEST_LIB NAMES armpl_lp64 armpl_lp64_mp armpl)

  find_package(nvpl QUIET)

  if(PkgConfig_FOUND)
    pkg_check_modules(_CRAY_LIBSCI_PKG QUIET sci)
  endif()

  # Do not guess a vendor whose BLA_VENDOR support requires a newer CMake than the one running:
  # fall back to the next candidate (eventually Generic) with a warning instead of failing later.
  set(_vpl_aocl_usable ${_AOCL_UTILS_LIB})
  if(_AOCL_UTILS_LIB AND CMAKE_VERSION VERSION_LESS 3.27)
    set(_vpl_aocl_usable FALSE)
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(WARNING "VendorPerfLibs: aoclutils library file found but AOCL support requires CMake 3.27 or newer (current: ${CMAKE_VERSION}). Ignoring AOCL.")
    endif()
  endif()
  set(_vpl_nvpl_usable ${nvpl_FOUND})
  if(nvpl_FOUND AND CMAKE_VERSION VERSION_LESS 4.1)
    set(_vpl_nvpl_usable FALSE)
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(WARNING "VendorPerfLibs: nvpl package found but NVPL support requires CMake 4.1 or newer (current: ${CMAKE_VERSION}). Ignoring NVPL.")
    endif()
  endif()

  if(_MKL_CORE_TEST_LIB)
    set(VPL_ID_GUESS "IntelMKL")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "Found mkl_core library file. Guessed VPL_ID 'IntelMKL'.")
    endif()
  elseif(_vpl_aocl_usable)
    set(VPL_ID_GUESS "AOCL")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "Found aoclutils library file. Guessed VPL_ID 'AOCL'.")
    endif()
  elseif(_ARMPL_TEST_LIB)
    set(VPL_ID_GUESS "ARMPL")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "Found armpl library file. Guessed VPL_ID 'ARMPL'.")
    endif()
  elseif(_vpl_nvpl_usable)
    set(VPL_ID_GUESS "NVPL")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "Found nvpl package. Guessed VPL_ID 'NVPL'.")
    endif()
  elseif(_CRAY_LIBSCI_PKG_FOUND)
    set(VPL_ID_GUESS "CrayLibSci")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "Found Cray libsci pkg-config package. Guessed VPL_ID 'CrayLibSci'.")
    endif()
  else()
    set(VPL_ID_GUESS "Generic")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "Signs of vendor performance libraries not found. Guessed VPL_ID 'Generic'.")
    endif()
  endif()
  unset(_vpl_aocl_usable)
  unset(_vpl_nvpl_usable)
endmacro()

if(NOT VPL_ID)
  if(BLA_VENDOR)
    # If BLA_VENDOR was set, interpret user intention as opting out of auto-detection by VPL.
    set(VPL_ID_GUESS "Generic")
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(STATUS "BLA_VENDOR has been set to '${BLA_VENDOR}'. Guessed VPL_ID 'Generic'.")
    endif()
  else()
    _vpl_speculate_vendor()
  endif()
else()
  _vpl_check_ID("VPL_ID" "${VPL_ID}")
  set(VPL_ID_GUESS "${VPL_ID}")
endif()

set(VPL_ID "${VPL_ID_GUESS}" CACHE STRING "Vendor Performance Library ID (IntelMKL, AOCL, ARMPL, NVPL, CrayLibSci, Generic)")
set_property(CACHE VPL_ID PROPERTY STRINGS ${_VPL_VALID_IDS})

set(_find_package_args)
if(VendorPerfLibs_FIND_QUIETLY)
  list(APPEND _find_package_args QUIET)
endif()

function(_vpl_warn_MKL_Fortran_ABI_issue CORE_LIB_VAR)
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

macro(_vpl_find_core)
  set(VPL_CORE_FOUND TRUE)
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
        ${_VPL_MKL_INC_HINTS}
      PATH_SUFFIXES mkl
    )

    if(BLAS_FOUND AND VendorPerfLibs_INCLUDE_DIR)
      _vpl_warn_MKL_Fortran_ABI_issue(BLAS_LIBRARIES)
      # Prefer the BLAS::BLAS target: unlike BLAS_LIBRARIES it also carries BLAS_LINKER_FLAGS.
      if(TARGET BLAS::BLAS)
        set(VPL_CORE_LIBRARIES BLAS::BLAS)
      else()
        set(VPL_CORE_LIBRARIES ${BLAS_LIBRARIES})
      endif()
    else()
      set(VPL_CORE_FOUND FALSE)
      if(NOT VendorPerfLibs_FIND_QUIETLY)
        set(_vpl_warning_core_mkl)
        if(NOT VendorPerfLibs_INCLUDE_DIR)
          list(APPEND _vpl_warning_core_mkl "Intel header file mkl.h not found.\n")
        endif()
        if(NOT BLAS_FOUND)
          # exactly the same check as in _vpl_speculate_vendor in case it was not called.
          find_library(_MKL_CORE_TEST_LIB NAMES mkl_core
            HINTS
              ${_VPL_MKL_LIB_HINTS}
          )
          if(_MKL_CORE_TEST_LIB)
            list(APPEND _vpl_warning_core_mkl "Intel mkl_core library file was found. FindBLAS failed most likely when testing linking MKL for complex reasons.\n")
            if(VPL_OMP)
              list(APPEND _vpl_warning_core_mkl "Threaded MKL requested but linking may fail due to incompatible compiler OpenMP runtime. To request sequential MKL, set VPL_OMP=OFF.\n")
            endif()
          else()
            list(APPEND _vpl_warning_core_mkl "Intel mkl_core library file not found. Please set the MKL root directory via a CMake variable CMAKE_PREFIX_PATH or VendorPerfLibs_ROOT.\n")
          endif()
        endif()
        message(WARNING ${_vpl_warning_core_mkl} "If you'd like to fully opt out Intel MKL, set VPL_ID=Generic.")
      endif()
    endif()
  elseif(VPL_ID STREQUAL "AOCL")
    if(CMAKE_VERSION VERSION_LESS 3.27)
      message(FATAL_ERROR "VendorPerfLibs: AOCL support requires CMake 3.27 or newer (current: ${CMAKE_VERSION})")
    endif()
    list(APPEND _vpl_required_vars VendorPerfLibs_INCLUDE_DIR)

    set(_vpl_warning_core_aocl)
    find_path(VendorPerfLibs_INCLUDE_DIR NAMES blis/blis.h blis.h)
    if(NOT VendorPerfLibs_INCLUDE_DIR)
      set(VPL_CORE_FOUND FALSE)
      list(APPEND _vpl_warning_core_aocl "AOCL header file blis.h not found.\n")
    endif()

    if(NOT VendorPerfLibs_FIND_QUIETLY AND NOT VendorPerfLibs_INCLUDE_DIR)
      message(WARNING ${_vpl_warning_core_aocl} "If you'd like to fully opt out AOCL, set VPL_ID=Generic.")
    endif()
  elseif(VPL_ID STREQUAL "ARMPL")
    list(APPEND _vpl_required_vars VendorPerfLibs_INCLUDE_DIR VPL_CORE_LIBRARIES)

    if(NOT VPL_OMP)
      set(BLA_VENDOR "Arm")
    else()
      set(BLA_VENDOR "Arm_mp")
    endif()
    find_package(BLAS ${_find_package_args})

    find_path(VendorPerfLibs_INCLUDE_DIR NAMES armpl.h)

    if(BLAS_FOUND AND VendorPerfLibs_INCLUDE_DIR)
      # Prefer the BLAS::BLAS target: unlike BLAS_LIBRARIES it also carries BLAS_LINKER_FLAGS.
      if(TARGET BLAS::BLAS)
        set(VPL_CORE_LIBRARIES BLAS::BLAS)
      else()
        set(VPL_CORE_LIBRARIES ${BLAS_LIBRARIES})
      endif()
    else()
      set(VPL_CORE_FOUND FALSE)
      if(NOT VendorPerfLibs_FIND_QUIETLY)
        set(_vpl_warning_core_armpl)
        if(NOT VendorPerfLibs_INCLUDE_DIR)
          list(APPEND _vpl_warning_core_armpl "ArmPL header file armpl.h not found.\n")
        endif()
        if(NOT BLAS_FOUND)
          list(APPEND _vpl_warning_core_armpl "ArmPL library not found.\n")
        endif()
        message(WARNING ${_vpl_warning_core_armpl} "If you'd like to fully opt out ArmPL, set VPL_ID=Generic.")
      endif()
    endif()
  elseif(VPL_ID STREQUAL "NVPL")
    if(CMAKE_VERSION VERSION_LESS 4.1)
      message(FATAL_ERROR "VendorPerfLibs: NVPL support requires CMake 4.1 or newer (current: ${CMAKE_VERSION})")
    endif()
    find_package(nvpl QUIET)
    if(NOT nvpl_FOUND)
      set(VPL_CORE_FOUND FALSE)
      if(NOT VendorPerfLibs_FIND_QUIETLY)
        message(WARNING "NVPL not found. If you'd like to fully opt out NVPL, set VPL_ID=Generic.")
      endif()
    endif()
  elseif(VPL_ID STREQUAL "CrayLibSci")
    list(APPEND _vpl_required_vars VendorPerfLibs_INCLUDE_DIR)

    if(VPL_OMP)
      set(_cray_libsci_name "sci_mp")
    else()
      set(_cray_libsci_name "sci")
    endif()

    set(BLA_PREFER_PKGCONFIG ON)
    set(BLA_PKGCONFIG_BLAS ${_cray_libsci_name})
    find_package(BLAS ${_find_package_args})
    if(BLAS_FOUND)
      set(BLAS_FIND_QUIETLY ${VendorPerfLibs_FIND_QUIETLY})
      find_package_message(BLAS "Found BLAS: ${BLAS_LIBRARIES}" "[${BLAS_LIBRARIES}]")
    endif()

    find_path(VendorPerfLibs_INCLUDE_DIR
      NAMES cblas.h lapacke.h
      HINTS
        ${PKGC_BLAS_INCLUDE_DIRS}
        ${BLAS_INCLUDE_DIRS}
      PATH_SUFFIXES include
    )

    if(NOT (BLAS_FOUND AND VendorPerfLibs_INCLUDE_DIR))
      set(VPL_CORE_FOUND FALSE)
      if(NOT VendorPerfLibs_FIND_QUIETLY)
        set(_vpl_warning_core_cray)
        if(NOT VendorPerfLibs_INCLUDE_DIR)
          list(APPEND _vpl_warning_core_cray "Cray libsci header file cblas.h not found.\n")
        endif()
        if(NOT BLAS_FOUND)
          list(APPEND _vpl_warning_core_cray "Cray libsci BLAS library not found.\n")
        endif()
        message(WARNING ${_vpl_warning_core_cray} "If you'd like to fully opt out Cray libsci, set VPL_ID=Generic.")
      endif()
    endif()
    unset(_cray_libsci_name)
  endif()

  if(VPL_CORE_FOUND)
    list(APPEND _vpl_lib_found_ids "core(${VPL_ID})")
  endif()
endmacro()

macro(_vpl_find_lapack)
  if(VPL_ID STREQUAL "IntelMKL")
    if(NOT VPL_OMP)
      set(BLA_VENDOR "Intel10_64lp_seq")
    else()
      set(BLA_VENDOR "Intel10_64lp")
    endif()
    find_package(LAPACK ${_find_package_args})
    if(LAPACK_FOUND)
      _vpl_warn_MKL_Fortran_ABI_issue(LAPACK_LIBRARIES)
    endif()
  elseif(VPL_ID STREQUAL "AOCL")
    if(NOT VPL_OMP)
      set(BLA_VENDOR "AOCL")
    else()
      set(BLA_VENDOR "AOCL_mt")
    endif()
    find_package(LAPACK ${_find_package_args})
  elseif(VPL_ID STREQUAL "ARMPL")
    if(NOT VPL_OMP)
      set(BLA_VENDOR "Arm")
    else()
      set(BLA_VENDOR "Arm_mp")
    endif()
    find_package(LAPACK ${_find_package_args})
  elseif(VPL_ID STREQUAL "NVPL")
    if(NOT VPL_OMP)
      set(BLA_THREAD "SEQ")
    else()
      set(BLA_THREAD "OMP")
    endif()
    set(BLA_VENDOR "NVPL")
    find_package(LAPACK ${_find_package_args})
  elseif(VPL_ID STREQUAL "CrayLibSci")
    if(VPL_OMP)
      set(_cray_libsci_name "sci_mp")
    else()
      set(_cray_libsci_name "sci")
    endif()
    set(BLA_PKGCONFIG_LAPACK ${_cray_libsci_name})
    find_package(LAPACK ${_find_package_args})
    if(LAPACK_FOUND)
      set(LAPACK_FIND_QUIETLY ${VendorPerfLibs_FIND_QUIETLY})
      find_package_message(LAPACK "Found LAPACK: ${LAPACK_LIBRARIES}" "[${LAPACK_LIBRARIES}]")
    endif()
    unset(_cray_libsci_name)
  else()
    find_package(LAPACK ${_find_package_args})
  endif()

  if(LAPACK_FOUND)
    # Prefer the LAPACK::LAPACK target: unlike LAPACK_LIBRARIES it also carries LAPACK_LINKER_FLAGS.
    if(TARGET LAPACK::LAPACK)
      set(VPL_LAPACK_LINK LAPACK::LAPACK)
    else()
      set(VPL_LAPACK_LINK ${LAPACK_LIBRARIES})
    endif()
    list(APPEND _vpl_lib_found_ids "lapack(${VPL_ID})")
    set(VendorPerfLibs_lapack_FOUND TRUE)
  else()
    set(VendorPerfLibs_lapack_FOUND FALSE)
    if(NOT VendorPerfLibs_FIND_QUIETLY)
      message(WARNING "LAPACK for VPL_ID '${VPL_ID}' with OpenMP ${VPL_OMP}, not found")
    endif()
  endif()
endmacro()

macro(_vpl_find_fft)
  # Default to VPL_ID without caching, so that changing VPL_ID later is not defeated by a stale cached value.
  if(NOT VPL_fft_ID)
    set(VPL_fft_ID "${VPL_ID}")
    list(APPEND _vpl_defaulted_ids VPL_fft_ID)
  endif()
  _vpl_check_ID("VPL_fft_ID" "${VPL_fft_ID}")

  set(VPL_FFT_FOUND TRUE)
  if(VPL_fft_ID STREQUAL "IntelMKL")
    if(NOT VPL_ID STREQUAL "IntelMKL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_fft_ID is IntelMKL but VPL_ID is not IntelMKL. Unsupported.")
    endif()
    find_path(VendorPerfLibs_FFTW_INCLUDE_DIR NAMES fftw3.f03
      HINTS
        ${_VPL_MKL_INC_HINTS}
      PATH_SUFFIXES fftw mkl/fftw
    )
    set(VPL_FFT_LIBRARIES ${VPL_CORE_LIBRARIES})
    if(NOT VendorPerfLibs_FFTW_INCLUDE_DIR OR NOT VPL_FFT_LIBRARIES)
      set(VPL_FFT_FOUND FALSE)
    endif()
  elseif(VPL_fft_ID STREQUAL "ARMPL")
    if(NOT VPL_ID STREQUAL "ARMPL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_fft_ID is ARMPL but VPL_ID is not ARMPL. Unsupported.")
    endif()
    find_path(VendorPerfLibs_FFTW_INCLUDE_DIR NAMES fftw3.f03 fftw3.h)
    set(VPL_FFT_LIBRARIES ${VPL_CORE_LIBRARIES})
    if(NOT VendorPerfLibs_FFTW_INCLUDE_DIR OR NOT VPL_FFT_LIBRARIES)
      set(VPL_FFT_FOUND FALSE)
    endif()
  elseif(VPL_fft_ID STREQUAL "NVPL")
    if(NOT VPL_ID STREQUAL "NVPL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_fft_ID is NVPL but VPL_ID is not NVPL. Unsupported.")
    endif()
    find_package(nvpl QUIET)
    if(TARGET nvpl::fftw)
      set(VPL_FFT_LIBRARIES nvpl::fftw)
    else()
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

macro(_vpl_find_vml)
  # Default to VPL_ID without caching, so that changing VPL_ID later is not defeated by a stale cached value.
  if(NOT VPL_vml_ID)
    set(VPL_vml_ID "${VPL_ID}")
    list(APPEND _vpl_defaulted_ids VPL_vml_ID)
  endif()
  _vpl_check_ID("VPL_vml_ID" "${VPL_vml_ID}")

  set(VPL_VML_FOUND TRUE)
  if(VPL_vml_ID STREQUAL "IntelMKL")
    if(NOT VPL_ID STREQUAL "IntelMKL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_vml_ID is IntelMKL but VPL_ID is not IntelMKL. Unsupported.")
    endif()
    # VML is part of MKL core. No additional libraries needed beyond core MKL libraries.
    set(VPL_VML_LIBRARIES ${VPL_CORE_LIBRARIES})
    if(NOT VPL_VML_LIBRARIES)
      set(VPL_VML_FOUND FALSE)
    endif()
  elseif(VPL_vml_ID STREQUAL "AOCL")
    if(NOT VPL_ID STREQUAL "AOCL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_vml_ID is AOCL but VPL_ID is not AOCL. Unsupported.")
    endif()
    find_library(VPL_VML_LIBRARIES NAMES alm)
    if(NOT VPL_VML_LIBRARIES)
      set(VPL_VML_FOUND FALSE)
    endif()
  elseif(VPL_vml_ID STREQUAL "ARMPL")
    if(NOT VPL_ID STREQUAL "ARMPL")
      message(FATAL_ERROR "VendorPerfLibs: VPL_vml_ID is ARMPL but VPL_ID is not ARMPL. Unsupported.")
    endif()
    find_library(VPL_VML_LIBRARIES NAMES amath)
    if(NOT VPL_VML_LIBRARIES)
      set(VPL_VML_FOUND FALSE)
    endif()
  else()
    # VML for NVPL, CrayLibSci, Generic is not currently supported (e.g. no direct open source drop-in)
    set(VPL_VML_FOUND FALSE)
  endif()

  if(VPL_VML_FOUND)
    list(APPEND _vpl_lib_found_ids "vml(${VPL_vml_ID})")
    set(VendorPerfLibs_vml_FOUND TRUE)
  else()
    set(VendorPerfLibs_vml_FOUND FALSE)
    # if VML is not supported by the vendor library or VPL_vml_ID=Generic, there is no need to issue the following warning.
    if(NOT VendorPerfLibs_FIND_QUIETLY AND DEFINED VPL_VML_LIBRARIES)
      message(WARNING "VML for VPL_vml_ID '${VPL_vml_ID}' not found")
    endif()
  endif()
endmacro()

set(_vpl_defaulted_ids)
unset(VendorPerfLibs_lapack_FOUND)
unset(VPL_CORE_LIBRARIES)
unset(VPL_LAPACK_LINK)
set(_vpl_lib_found_ids)
set(_vpl_required_vars _vpl_lib_found_ids)
_vpl_find_core()

if(NOT VendorPerfLibs_FIND_COMPONENTS OR "lapack" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
  _vpl_find_lapack()
endif()

if(NOT VendorPerfLibs_FIND_COMPONENTS OR "fft" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
  _vpl_find_fft()
endif()

if(NOT VendorPerfLibs_FIND_COMPONENTS OR "vml" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
  _vpl_find_vml()
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
        INTERFACE_LINK_LIBRARIES "${VPL_LAPACK_LINK}"
      )
      add_library(VPL::lapack ALIAS vpl_lapack)
    endif()
  endif()

  # Create vpl_fft
  if(NOT VendorPerfLibs_FIND_COMPONENTS OR "fft" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
    if(NOT TARGET vpl_fft)
      # Unquoted expansion drops empty variables (e.g. no vendor include dir for Generic).
      set(_vpl_fft_include_dirs ${VendorPerfLibs_INCLUDE_DIR} ${VendorPerfLibs_FFTW_INCLUDE_DIR})
      list(REMOVE_DUPLICATES _vpl_fft_include_dirs)
      add_library(vpl_fft INTERFACE IMPORTED)
      set_target_properties(vpl_fft PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_vpl_fft_include_dirs}"
        INTERFACE_LINK_LIBRARIES "${VPL_FFT_LIBRARIES}"
      )
      add_library(VPL::fft ALIAS vpl_fft)
      unset(_vpl_fft_include_dirs)
    endif()
  endif()

  # Create vpl_vml
  if(NOT VendorPerfLibs_FIND_COMPONENTS OR "vml" IN_LIST VendorPerfLibs_FIND_COMPONENTS)
    if(NOT TARGET vpl_vml)
      add_library(vpl_vml INTERFACE IMPORTED)
      set_target_properties(vpl_vml PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${VendorPerfLibs_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${VPL_VML_LIBRARIES}"
      )
      add_library(VPL::vml ALIAS vpl_vml)
    endif()
  endif()
endif()

# If the BLAS/LAPACK libraries were found, keep the BLA_* variables as set above so that later find_package(BLAS)
# and find_package(LAPACK) calls outside this module follow the same discovery.
# Otherwise restore the caller's original values.
# The same applies to the Cray prefixes appended to CMAKE_PREFIX_PATH.
# The lapack result is authoritative when the lapack search ran; otherwise (e.g. only fft or vml requested)
# use the result of the core (BLAS) search.
if(DEFINED VendorPerfLibs_lapack_FOUND)
  set(_vpl_blas_lapack_ok "${VendorPerfLibs_lapack_FOUND}")
else()
  set(_vpl_blas_lapack_ok "${VPL_CORE_FOUND}")
endif()
if(NOT _vpl_blas_lapack_ok)
  _vpl_restore_BLA_vars()
  set(CMAKE_PREFIX_PATH "${_vpl_saved_CMAKE_PREFIX_PATH}")
endif()
unset(_vpl_saved_CMAKE_PREFIX_PATH)
unset(_vpl_blas_lapack_ok)
_vpl_clear_BLA_vars()

cmake_pop_check_state()

# Scope cleanup
unset(_VPL_BLA_VARS)
unset(VPL_CORE_LIBRARIES)
unset(VPL_LAPACK_LINK)
unset(_VPL_MKL_LIB_HINTS)
unset(_VPL_MKL_INC_HINTS)
foreach(_vpl_id_var IN LISTS _vpl_defaulted_ids)
  unset(${_vpl_id_var})
endforeach()
unset(_vpl_id_var)
unset(_vpl_defaulted_ids)
unset(_find_package_args)
unset(_vpl_required_vars)
unset(_vpl_lib_found_ids)
unset(_VPL_VALID_IDS)
