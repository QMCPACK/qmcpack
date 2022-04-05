#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers
#//
#// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
#//
#// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
#//////////////////////////////////////////////////////////////////////////////////////

unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
set(SRC_DIR tools)
set(UTEST_EXE test_${SRC_DIR})
set(UTEST_NAME deterministic-unit_test_${SRC_DIR})
set(UTEST_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(SRCS test_qmcfstool.cpp)

execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple_Sk.dat" ${UTEST_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/simple_input.xml" ${UTEST_DIR})

add_executable(${UTEST_EXE} ${SRCS})
target_link_libraries(${UTEST_EXE} catch_main fstool)
if(USE_OBJECT_TARGET)
  target_link_libraries(${UTEST_EXE} qmcparticle qmcutil)
endif()

add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${UTEST_DIR})

# Minimal test for qmc-check-affinity. Check for "OpenMP" in output, no errorcodes.
if(HAVE_MPI)
  add_test(NAME check-affinity_runs COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                                            $<TARGET_FILE:qmc-check-affinity>)
else()
  add_test(NAME check-affinity_runs COMMAND $<TARGET_FILE:qmc-check-affinity>)
endif()

if(QMC_OMP)
  set_tests_properties(
    check-affinity_runs
    PROPERTIES TIMEOUT
               120
               LABELS
               "deterministic"
               PASS_REGULAR_EXPRESSION
               "OpenMP"
               PROCESSORS
               1
               ENVIRONMENT
               OMP_NUM_THREADS=3
               PROCESSOR_AFFINITY
               TRUE)
else()
  message(STATUS "Check-affinity-runs without OpenMP")
  set_tests_properties(check-affinity_runs PROPERTIES TIMEOUT 120 LABELS "deterministic" PROCESSORS 1)
endif()

# Basic test for qmc-get-supercell. Cubic primitive vectors, request what should optimally be a 2x2x2 tiling
if(HAVE_MPI)
  add_test(NAME get-supercell_verbose
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:qmc-get-supercell>
                   --ptvs 1 0 0 0 1 0 0 0 1 --target 8 --maxentry 4 --verbose)
else()
  add_test(NAME get-supercell_verbose COMMAND $<TARGET_FILE:qmc-get-supercell> --ptvs 1 0 0 0 1 0 0 0 1 --target 8
                                              --maxentry 4 --verbose)
endif()
set_tests_properties(
  get-supercell_verbose
  PROPERTIES TIMEOUT
             120
             LABELS
             "deterministic"
             PASS_REGULAR_EXPRESSION
             "radius = 1"
             PROCESSORS
             1
             ENVIRONMENT
             OMP_NUM_THREADS=3
             PROCESSOR_AFFINITY
             TRUE)
if(HAVE_MPI)
  add_test(NAME get-supercell_terse COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                                            $<TARGET_FILE:qmc-get-supercell> --ptvs 1 0 0 0 1 0 0 0 1 --target 8)
else()
  add_test(NAME get-supercell_terse COMMAND $<TARGET_FILE:qmc-get-supercell> --ptvs 1 0 0 0 1 0 0 0 1 --target 8)
endif()
set_tests_properties(
  get-supercell_terse
  PROPERTIES TIMEOUT
             120
             LABELS
             "deterministic"
             PASS_REGULAR_EXPRESSION
             "1   2   0   0   0   2   0   0   0   2   2   0   0   0   2   0   0   0   2"
             PROCESSORS
             1
             ENVIRONMENT
             OMP_NUM_THREADS=3
             PROCESSOR_AFFINITY
             TRUE)

# Check qmc-extract-eshdf-kvectors
set(KTEST_HDF_INPUT1 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp/pwscf.pwscf.h5)
set(KTEST_HDF_INPUT2 ${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp/pwscf.pwscf.h5)
if(HAVE_MPI)
  add_test(NAME extract-eshdf-kvectors_1x1x1
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:qmc-extract-eshdf-kvectors> ${KTEST_HDF_INPUT1})
else()
  add_test(NAME extract-eshdf-kvectors_1x1x1 COMMAND $<TARGET_FILE:qmc-extract-eshdf-kvectors> ${KTEST_HDF_INPUT1})
endif()
set_tests_properties(
  extract-eshdf-kvectors_1x1x1
  PROPERTIES TIMEOUT
             120
             LABELS
             "deterministic"
             PASS_REGULAR_EXPRESSION
             "0   0   0"
             PROCESSORS
             1
             ENVIRONMENT
             OMP_NUM_THREADS=3
             PROCESSOR_AFFINITY
             TRUE)
if(HAVE_MPI)
  add_test(NAME extract-eshdf-kvectors_2x1x1
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                   $<TARGET_FILE:qmc-extract-eshdf-kvectors> ${KTEST_HDF_INPUT2})
else()
  add_test(NAME extract-eshdf-kvectors_2x1x1 COMMAND $<TARGET_FILE:qmc-extract-eshdf-kvectors> ${KTEST_HDF_INPUT2})
endif()
set_tests_properties(
  extract-eshdf-kvectors_2x1x1
  PROPERTIES TIMEOUT
             120
             LABELS
             "deterministic"
             PASS_REGULAR_EXPRESSION
             "0   0   0\n0.5   0   0"
             PROCESSORS
             1
             ENVIRONMENT
             OMP_NUM_THREADS=3
             PROCESSOR_AFFINITY
             TRUE)

# Minimal test for convert4qmc. Check for Usage on output if no arguments. Primarily tested through converter tests.
if(HAVE_MPI)
  add_test(NAME convert4qmc_usage COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                                          $<TARGET_FILE:convert4qmc>)
else()
  add_test(NAME convert4qmc_usage COMMAND $<TARGET_FILE:convert4qmc>)
endif()
set_tests_properties(
  convert4qmc_usage
  PROPERTIES TIMEOUT
             120
             LABELS
             "deterministic"
             PASS_REGULAR_EXPRESSION
             "Usage:"
             PROCESSORS
             1
             ENVIRONMENT
             OMP_NUM_THREADS=3
             PROCESSOR_AFFINITY
             TRUE)
