
cmake_minimum_required(VERSION 3.6.0)

enable_testing()
include(CTest)

SET(FAKE_MPI_DIR "..")
SET(MPI_WRAPPER_DIR "../../..")

set(TEST_SRCS
 test_c_mpi_init.c
 test_mpi_init.cpp
 test_wrapper_mpi_init.cpp
)

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}x.x")
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_include_directories(${TEST_EXE} PUBLIC "..")
  if (TEST_FILE MATCHES ".cpp$")
    target_include_directories(${TEST_EXE} PUBLIC ${MPI_WRAPPER_DIR})
  endif()
  add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
endforeach()


# Compile-time check for multiply defined symbols
add_library(library_check ${MPI_WRAPPER_DIR}/test/library_check.cpp)
add_executable(library_main ${MPI_WRAPPER_DIR}/test/library_main.cpp)
target_link_libraries(library_main library_check)
target_include_directories(library_check PUBLIC ${FAKE_MPI_DIR} ${MPI_WRAPPER_DIR})
target_include_directories(library_main PUBLIC ${FAKE_MPI_DIR} ${MPI_WRAPPER_DIR})

