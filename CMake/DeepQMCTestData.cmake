#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2026 QMCPACK developers.
#//////////////////////////////////////////////////////////////////////////////////////
# Local-only, trusted DeepQMC He fixture validation and test environment.
get_filename_component(QMC_DEEPQMC_TEST_DATA_DIR "${QMC_DEEPQMC_TEST_DATA_DIR}" ABSOLUTE)
set(_deepqmc_manifest "${qmcpack_SOURCE_DIR}/tests/molecules/He_deepqmc/SHA256SUMS")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${_deepqmc_manifest}"
  "${QMC_DEEPQMC_TEST_DATA_DIR}/chkpt-30000.pt"
  "${QMC_DEEPQMC_TEST_DATA_DIR}/config.yaml")
file(STRINGS "${_deepqmc_manifest}" _deepqmc_hashes)
foreach(_required IN ITEMS "chkpt-30000.pt" "config.yaml")
  if(NOT EXISTS "${QMC_DEEPQMC_TEST_DATA_DIR}/${_required}")
    message(FATAL_ERROR "Missing ${_required} in QMC_DEEPQMC_TEST_DATA_DIR=${QMC_DEEPQMC_TEST_DATA_DIR}; extract the trusted He fixture locally and configure again")
  endif()
  set(_found_hash "")
  foreach(_entry IN LISTS _deepqmc_hashes)
    if(_entry MATCHES "^([0-9a-f]+) [ *]${_required}$")
      set(_found_hash "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(NOT _found_hash)
    message(FATAL_ERROR "${_deepqmc_manifest} has no checksum entry for ${_required}")
  endif()
  file(SHA256 "${QMC_DEEPQMC_TEST_DATA_DIR}/${_required}" _actual_hash)
  if(NOT _actual_hash STREQUAL _found_hash)
    message(FATAL_ERROR "Checksum mismatch for ${QMC_DEEPQMC_TEST_DATA_DIR}/${_required}; use the trusted He fixture matching ${_deepqmc_manifest}")
  endif()
endforeach()
set(QMC_DEEPQMC_HE_CHECKPOINT "${QMC_DEEPQMC_TEST_DATA_DIR}/chkpt-30000.pt")
set(QMC_DEEPQMC_HYDRA_CONFIG "${QMC_DEEPQMC_TEST_DATA_DIR}/config.yaml")
execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('purelib')); print(sysconfig.get_path('platlib'))"
  RESULT_VARIABLE _python_path_result OUTPUT_VARIABLE _python_paths OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT _python_path_result EQUAL 0)
  message(FATAL_ERROR "Unable to query Python site-package paths using ${Python3_EXECUTABLE}")
endif()
string(REPLACE "\n" ";" _python_paths "${_python_paths}")
list(GET _python_paths 0 QMC_DEEPQMC_PURELIB)
list(GET _python_paths 1 QMC_DEEPQMC_PLATLIB)
function(deepqmc_set_test_environment test_name)
  set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT
    "DEEPQMC_HE_CHECKPOINT=${QMC_DEEPQMC_HE_CHECKPOINT};DEEPQMC_HYDRA_CONFIG=${QMC_DEEPQMC_HYDRA_CONFIG};DEEPQMC_PYTHON_SITE_PACKAGES=${QMC_DEEPQMC_PURELIB};JAX_PLATFORMS=cpu;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
  set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT_MODIFICATION
    "PYTHONPATH=path_list_prepend:${QMC_DEEPQMC_PLATLIB};PYTHONPATH=path_list_prepend:${QMC_DEEPQMC_PURELIB}")
endfunction()
