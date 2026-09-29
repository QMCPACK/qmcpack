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

# FindFFTW.cmake
#
# Finds the FFTW library.
#
# This module supports the following components:
#   seq - sequential FFTW library (fftw3)
#   omp - OpenMP FFTW library (fftw3_omp)
#
# This module defines the following variables:
#   FFTW_FOUND        - True if FFTW was found.
#   FFTW_INCLUDE_DIR  - Include directory for FFTW.
#   FFTW_LIBRARIES    - All FFTW libraries found.
#
# This module creates the following imported targets:
#   FFTW::seq - The sequential FFTW library.
#   FFTW::omp - The OpenMP FFTW library.

include(FindPackageHandleStandardArgs)

# Find include directory
find_path(FFTW_INCLUDE_DIR
  NAMES fftw3.h fftw3.f03
)

set(FFTW_LIBRARIES)
set(_FFTW_REQUIRED_VARS FFTW_INCLUDE_DIR)

# Find seq component
if(NOT FFTW_FIND_COMPONENTS OR "seq" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(FFTW_seq_LIBRARY NAMES fftw3)
  if(FFTW_seq_LIBRARY)
    set(FFTW_seq_FOUND TRUE)
    list(PREPEND FFTW_LIBRARIES ${FFTW_seq_LIBRARY})
  else()
    set(FFTW_seq_FOUND FALSE)
  endif()
  list(PREPEND _FFTW_REQUIRED_VARS FFTW_seq_LIBRARY)
endif()

# Find omp component
if(NOT FFTW_FIND_COMPONENTS OR "omp" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(FFTW_omp_LIBRARY NAMES fftw3_omp)
  if(FFTW_omp_LIBRARY)
    set(FFTW_omp_FOUND TRUE)
    list(PREPEND FFTW_LIBRARIES ${FFTW_omp_LIBRARY})
  else()
    set(FFTW_omp_FOUND FALSE)
  endif()
  list(PREPEND _FFTW_REQUIRED_VARS FFTW_omp_LIBRARY)
endif()

find_package_handle_standard_args(FFTW
  REQUIRED_VARS ${_FFTW_REQUIRED_VARS}
  HANDLE_COMPONENTS
)
