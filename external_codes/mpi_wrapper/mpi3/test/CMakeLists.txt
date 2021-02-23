# Use cmake/ctest for building and testing boost mpi3 tests

# Not sure of the actual required minimum.
cmake_minimum_required(VERSION 3.6.0)
project(boostmpi3_tests)

enable_testing()
include(CTest)

if (NOT ($ENV{CRAYPE_VERSION} MATCHES "."))
  set(CMAKE_C_COMPILER mpicc)
  set(CMAKE_CXX_COMPILER mpic++)
endif()

# Find the correct parallel program launcher
if ($ENV{CRAYPE_VERSION} MATCHES ".")
  find_program(MPIEXEC aprun)

  # Using Slurm
  if (NOT MPIEXEC)
     find_program(MPIEXEC srun)
  endif()
endif()
# Everywhere else likely uses mpirun
if (NOT MPIEXEC)
  set(MPIEXEC mpirun)
endif()


# Set BOOST_ROOT if boost is in a non-standard location
find_package(Boost REQUIRED)
if(Boost_FOUND)
  set(HAVE_LIBBOOST 1)
  include_directories(${Boost_INCLUDE_DIRS})
  message(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")
endif()

# This list is only a subset of files.  It might be useful to divide some of the tests into
#  different categories (MPI-1, MPI-2, etc., or something else)
set(TEST_SRCS
communicator_abort.cpp
communicator_all_gather.cpp
communicator_all_gatherv.cpp
communicator_all_gatherv_output_iterator.cpp
communicator_all_reduce.cpp
communicator_barrier.cpp
communicator_broadcast.cpp
communicator_divide.cpp
gather2.cpp
communicator_ibroadcast.cpp
communicator_igather.cpp
communicator_list.cpp
communicator_operator.cpp
communicator_reduce.cpp
communicator_reduce_in_place.cpp
communicator_scatter.cpp
communicator_set_error_handler.cpp
communicator_split.cpp
window_put_get.cpp
shared_mutex.cpp
)

# List of tests that need to link with boost_serialization
SET(NEED_BOOST_SERIALIZATION_SRCS
   communicator_send.cpp
   communicator_send_class_nonintrusive.cpp
   communicator_send_receive.cpp
   communicator_send_class.cpp
   communicator_ireceive.cpp
   communicator_gather.cpp
#   process.cpp
#   process_vector.cpp
)

# Build tests that need boost serialization
set(TEST_BOOST_SERIALIZATION TRUE)
if (TEST_BOOST_SERIALIZATION)
  SET(TEST_SRCS ${TEST_SRCS} ${NEED_BOOST_SERIALIZATION_SRCS})
endif()

set(CMAKE_CXX_STANDARD 14)

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}.x")
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_compile_options(${TEST_EXE} PUBLIC -Wall -Wextra -Wpedantic -Werror -Wfatal-errors)
  target_include_directories(${TEST_EXE} PUBLIC "../../../..")

  list(FIND NEED_BOOST_SERIALIZATION_SRCS ${TEST_FILE} NEED_BOOST_SERIALIZATION)
  if (NOT (${NEED_BOOST_SERIALIZATION} EQUAL -1))
    target_compile_options(${TEST_EXE} PUBLIC -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY)
#    target_link_libraries(${TEST_EXE} boost_serialization)
  endif()

  target_link_libraries(${TEST_EXE} pthread)

  set(NPROC 3)
  if ((TEST_FILE STREQUAL "communicator_divide.cpp"  ) OR
      (TEST_FILE STREQUAL "communicator_operator.cpp") OR
      (TEST_FILE STREQUAL "communicator_scatter.cpp" ))
    set(NPROC 6)
  endif()
  if(DEFINED ENV{VALGRIND_EXE})
    add_test(NAME ${TEST_EXE} COMMAND ${MPIEXEC} $ENV{MPI_OVERSUBSCRIBE} -n ${NPROC} valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all --error-exitcode=1 ./${TEST_EXE})
  else()
    add_test(NAME ${TEST_EXE} COMMAND ${MPIEXEC} $ENV{MPI_OVERSUBSCRIBE} -n ${NPROC} ./${TEST_EXE})
  endif()
endforeach()

set_tests_properties(communicator_abort.cpp.x PROPERTIES WILL_FAIL true)

# Compile-time check for multiply defined symbols
add_library(library_check library_check.cpp)
add_executable(library_main.x library_main.cpp)
target_link_libraries(library_main.x library_check)

