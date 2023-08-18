# Use cmake/ctest for building and testing boost mpi3 tests
cmake_minimum_required(VERSION 3.18.4) # 3.10 for FindMPI, 3.12 for MPI_ROOT
project(boostmpi3_tests VERSION 0.1 LANGUAGES CXX)

if(EXISTS "CMakeLists.txt")
	message(FATAL_ERROR "You probably don't want to run CMake in a directory with CMakeLists.txt")
endif()

enable_testing()
include(CTest)
include(FetchContent)

FetchContent_Declare(multi GIT_REPOSITORY https://gitlab.com/correaa/boost-multi.git)  # needs git installed
# FetchContent_Declare(multi URL https://gitlab.com/correaa/boost-multi/-/archive/master/boost-multi-master.zip DOWNLOAD_EXTRACT_TIMESTAMP ON)
FetchContent_MakeAvailable(multi)

find_package(MPI REQUIRED)  # install MPI, e.g. `sudo dnf install mpich-devel` or `module load mpi`

include_directories(SYSTEM ${MPI_CXX_INCLUDE_DIRS})
link_libraries(${MPI_CXX_LIBRARIES})

if (NOT ($ENV{CRAYPE_VERSION} MATCHES "."))
#  set(CMAKE_C_COMPILER mpicc)
#  set(CMAKE_CXX_COMPILER mpic++)
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
#  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
#  message(STATUS "Setting Boost_INCLUDE_DIRS=${Boost_INCLUDE_DIRS}")

  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  message(STATUS "Setting Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
endif()

if(ENABLE_CUDA)
	enable_language(CUDA)
	find_package(CUDA 11.3.0 REQUIRED)
#	set(CXX_COMPILER ${CMAKE_CUDA_COMPILER})
#	set(CXX_FLAGS ${CMAKE_CUDA_FLAGS})
#	string(APPEND CMAKE_CUDA_FLAGS " --forward-unknown-to-host-linker") # ${ARCH_FLAGS} -std=c++17 --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--diag_suppress=implicit_return_from_non_void_function\"")
#	add_compile_options(--forward-unknown-to-host-linker)
endif()

# This list is only a subset of files.  It might be useful to divide some of the tests into
#  different categories (MPI-1, MPI-2, etc., or something else)
set(TEST_SRCS
  all_reduce.cpp
	async_interaction.cpp
  broadcast.cpp
	cartesian.cpp
	communicator_abort.cpp
	communicator_all_gather.cpp
	communicator_all_gatherv.cpp
	communicator_all_gatherv_output_iterator.cpp
	communicator_barrier.cpp
	communicator_divide.cpp
	datatype.cpp
  datatype_struct_vector3.cpp
  deino_all_to_all.cpp
  deino_broadcast.cpp
  deino_op_create.cpp
	gather2.cpp
	group.cpp
	communicator_cctor.cpp
	communicator_ibroadcast.cpp
	communicator_igather.cpp
	communicator_grip_handle.cpp
	communicator_list.cpp
	communicator_main.cpp
	communicator_mutable.cpp
	communicator_operator.cpp
	communicator_ostream.cpp
	communicator_passby.cpp
	communicator_reduce.cpp
	communicator_reduce_in_place.cpp
	communicator_scatter.cpp
	communicator_set_error_handler.cpp
	communicator_split.cpp
	ibarrier.cpp
	environment_thread.cpp
	empty_main.cpp
	window_put_get.cpp
	reduce_maxloc.cpp
	ring.cpp
	shared_mutex.cpp
  send_complex.cpp
	simple_send_receive.cpp
	simple_broadcast.cpp
  spinor.cpp
  status.cpp
	type_commit.cpp
	type_size.cpp
	uniform_abort.cpp
)

# List of tests that need to link with boost_serialization
SET(NEED_BOOST_SERIALIZATION_SRCS
	communicator_send.cpp
	communicator_send_class_nonintrusive.cpp
	communicator_send_receive.cpp
	communicator_send_class.cpp
	communicator_ireceive.cpp
	communicator_gather.cpp
  process.cpp
#   process_vector.cpp
)

# Build tests that need boost serialization
set(TEST_BOOST_SERIALIZATION TRUE)
if (TEST_BOOST_SERIALIZATION)
	SET(TEST_SRCS ${TEST_SRCS} ${NEED_BOOST_SERIALIZATION_SRCS})
endif()

foreach(TEST_FILE ${TEST_SRCS})
  SET(TEST_EXE "${TEST_FILE}.x")
  add_executable(${TEST_EXE} ${TEST_FILE})
  target_compile_features(${TEST_EXE} PRIVATE cxx_std_17)

  if(ENABLE_CUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
  endif()

  target_compile_options(${TEST_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Werror -Wno-error=deprecated-declarations>) # -Wpedantic -Wunused -Wnon-virtual-dtor -Woverloaded-virtual -Wcast-qual -Wformat=2>) #-Wsign-conversion -Wshadow -Wconversion
  target_compile_options(${TEST_EXE} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Wall -Wextra>)

  target_compile_options(
    ${TEST_EXE}
    PRIVATE
      $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
        -Xcompiler=-Werror,-Wall,-Wextra,-Wcast-align,-Wcast-qual,-Wno-deprecated-declarations,-Wno-double-promotion,-Wduplicated-branches,-Wduplicated-cond,-Wformat-truncation,-Wformat=2,-Wlogical-op,-Wmisleading-indentation,-Wno-missing-include-dirs,-Wnon-virtual-dtor,-Wno-missing-declarations,-Wnon-virtual-dtor,-Wnull-dereference,-Woverloaded-virtual,-Wpointer-arith,-Wno-redundant-decls,-Wno-shadow,-Wno-switch-enum,-Wno-unknown-pragmas,-Wtrampolines,-Wuninitialized,-Wunused,-Wunused-but-set-variable,-Wunused-result,-Wno-zero-as-null-pointer-constant
        --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe=--display_error_number -Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier
      >
      $<$<AND:$<CXX_COMPILER_ID:GNU>,$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>,$<NOT:$<CUDA_COMPILER_ID:Clang>>>:
        -Wno-unknown-pragmas  
        # -Wextra-semi (gcc 8, not in gcc 7)
        -Wpedantic
        # -Wabi=13 -Wabi-tag (maybe important when linking with very old libraries)
        -Waddress
        # -Waddress-of-packed-member (gcc 11, not in gcc 8)
        # -Waggregate-return (disallow return classes or structs, seems a C-compatibility warning)
        -Waggressive-loop-optimizations
        -Wzero-as-null-pointer-constant
      >
      $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CUDA_COMPILER_ID:Clang>>:
        -Wno-unknown-pragmas
        -Wcast-align
        -Wcomplex-component-init
        -Wconstexpr-not-const
        -Wconsumed
        -Wconversion -Wconversion-null
        -Wmove
        -Wzero-as-null-pointer-constant -Wzero-length-array
        -Wno-c++98-compat-unnamed-type-template-args
        -Wno-ignored-qualifiers
        -Wno-range-loop-analysis # (this check is overzealous in clang 9)
      >
      $<$<CXX_COMPILER_ID:Intel>:  # also IntelLLVM, XL (ibm), XLClang (ibm)
        -Werror
        -Wall
        -Wextra
        -diag-disable=remark
        -wd161
        -Wabi
        -Warray-bounds
        -Wcast-qual
        -Wchar-subscripts
        -Wunused -Wunused-but-set-variable -Wunused-function -Wunused-parameter -Wunused-variable
        -Wwrite-strings
      >
      $<$<OR:$<CXX_COMPILER_ID:PGI>,$<CXX_COMPILER_ID:NVHPC>>:
        -Werror
        -Wall
        -Wcast-qual
        -Wformat=2
        -Wshadow
        -Wunused-parameter
      >
      $<$<CXX_COMPILER_ID:MSVC>:
        /WX
        /W4
      >
  )

  target_include_directories(${TEST_EXE} PUBLIC "../../../..")

  # TODO(correaa) simplify include directories
  target_include_directories(${TEST_EXE} PUBLIC "../include")
  target_include_directories(${TEST_EXE} PUBLIC "../include/mpi3")
  target_include_directories(${TEST_EXE} PUBLIC "../include/mpi3/dummy")

  #target_link_libraries(${TEST_EXE} PRIVATE mpi3)

  list(FIND NEED_BOOST_SERIALIZATION_SRCS ${TEST_FILE} NEED_BOOST_SERIALIZATION)
  if (NOT (${NEED_BOOST_SERIALIZATION} EQUAL -1))
    target_compile_options(${TEST_EXE} PUBLIC -D_MAKE_BOOST_SERIALIZATION_HEADER_ONLY)
#    target_link_libraries(${TEST_EXE} PRIVATE boost_serialization)
  endif()

  target_link_libraries(${TEST_EXE} PRIVATE pthread)
  target_link_libraries(${TEST_EXE} PRIVATE multi)

  set(NPROC 3)
  if (
	   (TEST_FILE STREQUAL "communicator_divide.cpp"     )
    OR (TEST_FILE STREQUAL "communicator_operator.cpp"   )
    OR (TEST_FILE STREQUAL "communicator_scatter.cpp"    )
    OR (TEST_FILE STREQUAL "cartesian.cpp"               )
	OR (TEST_FILE STREQUAL "simple_send_receive.cpp"     )
  )
    set(NPROC 6)
  endif()

  if (
       (TEST_FILE STREQUAL "uniform_abort.cpp"           )
    OR (TEST_FILE STREQUAL "communicator_grip_handle.cpp")
  )
    set(NPROC 4)
  endif()

  if (
       (TEST_FILE STREQUAL "cartesian.cpp"            )
  )
    set(NPROC 6)
  endif()

  if (
	   (TEST_FILE STREQUAL "reduce_maxloc.cpp"  )
  )
    set(NPROC 4)
  endif()

  if(DEFINED ENV{VALGRIND_EXE})
    add_test(NAME ${TEST_EXE} COMMAND ${MPIEXEC} $ENV{MPI_OVERSUBSCRIBE} -n ${NPROC} valgrind --leak-check=full --track-origins=yes --show-leak-kinds=all --error-exitcode=1 ./${TEST_EXE})
  else()
    add_test(NAME ${TEST_EXE} COMMAND ${MPIEXEC} $ENV{MPI_OVERSUBSCRIBE} -n ${NPROC} ./${TEST_EXE})
  endif()
  set_tests_properties(${TEST_EXE} PROPERTIES TIMEOUT 30)
endforeach()

set_tests_properties(communicator_abort.cpp.x PROPERTIES WILL_FAIL true)
set_tests_properties(uniform_abort.cpp.x      PROPERTIES WILL_FAIL true)

# Compile-time check for multiply defined symbols
add_library(library_check library_check.cpp)
target_compile_features(library_check PUBLIC cxx_std_17)

# TODO(correaa) simplify include directories
target_include_directories(library_check PUBLIC "../include")
target_include_directories(library_check PUBLIC "../include/mpi3")
target_include_directories(library_check PUBLIC "../include/mpi3/dummy")

add_executable(library_main.x library_main.cpp)
target_compile_features(library_main.x PRIVATE cxx_std_17)
target_link_libraries(library_main.x library_check)
