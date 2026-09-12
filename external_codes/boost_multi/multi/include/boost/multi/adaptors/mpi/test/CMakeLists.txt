cmake_minimum_required(VERSION 3.15)

if(ENABLE_CUDA)
	enable_language(CUDA)
endif()

enable_testing()
include(CTest)

find_package(Boost REQUIRED COMPONENTS NO_MODULE)

add_executable(mpi.cpp.x mpi.cpp)

if(ENABLE_CUDA)
    set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
endif()

if(APPLE)
    # https://apple.stackexchange.com/a/121010
    # add_custom_command (TARGET mpi.cpp.x POST_BUILD COMMAND codesign --force --deep -s ACTMPI ${CMAKE_CURRENT_BINARY_DIR}/mpi.cpp.x)
endif()

target_link_libraries(mpi.cpp.x PRIVATE multi)
target_link_libraries(mpi.cpp.x PRIVATE MPI::MPI_CXX Boost::boost)

if(APPLE)
    add_test(NAME mpi.cpp.x COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe --mca btl ^tcp ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi.cpp.x> ${MPIEXEC_POSTFLAGS})
else()
    add_test(NAME mpi.cpp.x COMMAND ${MPIEXEC_EXECUTABLE}                ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi.cpp.x> ${MPIEXEC_POSTFLAGS})
endif()