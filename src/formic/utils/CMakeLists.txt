# -std=c++11 is handled at the top level. Ye
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

REMOVE_DEFINITIONS(-DMPICH_SKIP_MPICXX)

# Really ugly hack to see if we're using MKL - replace it
#STRING(FIND "${BLAS_LIBRARIES}" "mkl" apos)
#MESSAGE("BLAS_LIBRARIES FIND = ${apos}")
#IF (NOT apos EQUAL -1)
#    MESSAGE("mkl found")
#    SET(MKL_FOUND 1)
#ENDIF()

# determine how to declare and how to mangle the blas/lapack functions
if(MKL_FOUND)
  set(LAPACK_DECLARE_INC formic/utils/declare_mkl.h)
  set(LAPACK_MANGLE_INC formic/utils/mangle_mkl.h)
else(MKL_FOUND)
  # use the fortran-C interface mangling determined by cmake
  # this should be correct so long as the same compilers were used to build blas/lapack
  # FortranCInterface_HEADER(${FORMIC_BINARY_DIR}/src/formic/utils/mangle_other.h MACRO_NAMESPACE "FC_")
  set(LAPACK_DECLARE_INC formic/utils/declare_other.h)
  set(LAPACK_MANGLE_INC formic/utils/mangle_other.h)
endif(MKL_FOUND)

# create header for declaration of blas/lapack functions
configure_file(declare.h.in ${FORMIC_BINARY_DIR}/src/formic/utils/declare.h)

# create header for mangling of blas/lapack function names
configure_file(mangle.h.in ${FORMIC_BINARY_DIR}/src/formic/utils/mangle.h)

# determine how to comment/uncomment the mpi interface header file
if(USE_MPI)
  set(MPI_COMMENT_I "")
  set(MPI_COMMENT_O "//")
else(USE_MPI)
  set(MPI_COMMENT_I "//")
  set(MPI_COMMENT_O "")
endif(USE_MPI)

# determine the level of mpi thread support to demand
if(USE_OPENMP)
  set(MPI_REQUESTED_THREADING "MPI_THREAD_SERIALIZED")
  #set(MPI_REQUESTED_THREADING "MPI_THREAD_FUNNELED")
  #set(MPI_REQUESTED_THREADING "MPI_THREAD_SINGLE")
else(USE_OPENMP)
  set(MPI_REQUESTED_THREADING "MPI_THREAD_SERIALIZED")
  #set(MPI_REQUESTED_THREADING "MPI_THREAD_SINGLE")
endif(USE_OPENMP)

# create header file for mpi functions
configure_file(mpi_interface.h.in ${FORMIC_BINARY_DIR}/src/formic/utils/mpi_interface.h)

# create cpp file for mpi functions
configure_file(mpi_interface.cpp.in ${FORMIC_BINARY_DIR}/src/formic/utils/mpi_interface.cpp)

# ensure we can use the eigen headers
include_directories(${EIGEN_INCLUDE_PATH})

add_library(formic_utils matrix.cpp
                         numeric.cpp
                         numeric_matrix.cpp
                         porabola_min_max.cpp
                         reusable_array.cpp
                         lapack_interface.cpp
                         random.cpp
                         timing.cpp
                         ${FORMIC_BINARY_DIR}/src/formic/utils/mpi_interface.cpp
                         lmyengine/davidson_solver.cpp
                         lmyengine/block_mat.cpp
                         lmyengine/block_alg.cpp
                         lmyengine/block_detail.cpp
                         lmyengine/energy_target.cpp
                         lmyengine/energy_target_accu.cpp
                         lmyengine/engine.cpp
                         lmyengine/engine_numeric.cpp
                         lmyengine/engine_timing.cpp
                         lmyengine/eom.cpp
                         lmyengine/matrix_builder.cpp
                         lmyengine/spam_solver.cpp
                         lmyengine/updater.cpp
                         lmyengine/var_dependencies.cpp)

#target_link_libraries(formic_utils ${FORMIC_LAPACK_BLAS_LIBS})

#set_property(TARGET formic_utils PROPERTY CXX_STANDARD 11)
#set_property(TARGET formic_utils PROPERTY CXX_STANDARD_REQUIRED ON)

#file(MAKE_DIRECTORY ${FORMIC_BINARY_DIR}/lib)
#file(MAKE_DIRECTORY ${FORMIC_BINARY_DIR}/include)
#file(MAKE_DIRECTORY ${FORMIC_BINARY_DIR}/include/formic)
#file(MAKE_DIRECTORY ${FORMIC_BINARY_DIR}/include/formic/utils)
#file(MAKE_DIRECTORY ${FORMIC_BINARY_DIR}/include/formic/utils/lmyengine)

#add_custom_command( TARGET formic_utils
#                    POST_BUILD
#                    COMMAND cp ${FORMIC_BINARY_DIR}/src/formic/utils/libformic_utils.* ${FORMIC_BINARY_DIR}/lib
#                    COMMAND cp ${FORMIC_BINARY_DIR}/src/formic/utils/*.h ${FORMIC_BINARY_DIR}/include/formic/utils
#                    COMMAND cp ${FORMIC_SOURCE_DIR}/src/formic/utils/*.h ${FORMIC_BINARY_DIR}/include/formic/utils
#                    COMMAND cp ${FORMIC_SOURCE_DIR}/src/formic/utils/lmyengine/*.h ${FORMIC_BINARY_DIR}/include/formic/utils/lmyengine
#                    COMMAND ${CMAKE_COMMAND} -E echo ""
#                    COMMAND ${CMAKE_COMMAND} -E echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#                    COMMAND ${CMAKE_COMMAND} -E echo "!!  The formic utilities library has been built and is located in lib.  !!"
#                    COMMAND ${CMAKE_COMMAND} -E echo "!!  The associated headers are located in include.                      !!"
#                    COMMAND ${CMAKE_COMMAND} -E echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
#                    COMMAND ${CMAKE_COMMAND} -E echo ""
#                    WORKING_DIRECTORY ${FORMIC_BINARY_DIR}
#                  )

#install(
#        TARGETS
#          formic_utils
#        RUNTIME DESTINATION
#          bin
#        LIBRARY DESTINATION
#          lib
#        ARCHIVE DESTINATION
#          lib
#       )
