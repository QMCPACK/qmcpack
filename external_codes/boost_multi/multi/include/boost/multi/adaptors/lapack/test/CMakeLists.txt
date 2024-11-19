cmake_minimum_required(VERSION 3.5)

if(ENABLE_CUDA)
	enable_language(CUDA)
endif()

# project(
#   boost-multi-adaptors-lapack-test
#   VERSION 0.1
#   LANGUAGES CXX
# )

# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Boost REQUIRED)
# add_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
# include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
# link_libraries(${Boost_LIBRARIES})
# link_directories(${Boost_LIBRARY_DIRS})

#find_package(LAPACK REQUIRED)

# find_package(BLAS REQUIRED)

# set(BLA_VENDOR OpenBLAS)
# set(BLA_VENDOR Intel10_64lp) find_package(BLAS) if(BLAS_FOUND) # in some systems with MKL, regular BLAS headers need to be found for it to work message("Multi/BLAS: MKL environment detected") add_definitions(-DRETURN_BY_STACK) else()
# message("Multi/BLAS: MKL environment not detected, looking for other BLAS") unset(BLA_VENDOR)
# find_package(BLAS REQUIRED)
# endif()

# find_path(
#   BLAS_INCLUDE_DIRS
#   cblas.h
#   /usr/include
#   /usr/local/include
#   $ENV{BLAS_HOME}/include
# )

# include_directories(../../../../../include)

# link_libraries(${BLAS_LIBRARIES})
#link_libraries(-llapacke)

# include_directories(${TEST_EXE} PRIVATE ${BLAS_INCLUDE_DIRS})

# add_compile_options(
#   -Werror
#   -Wall
#   -Wextra
#   -fno-common
#   -Wfatal-errors
#   # $<$<CXX_COMPILER_ID:GNU>: -Wpedantic -Wformat-truncation -fstack-usage >#-Wconversion
#   # $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
#   # -Wpedantic -Wmove > $<$<CXX_COMPILER_ID:Intel>: -wd161 -diag-disable=remark -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses
#   # -Wpointer-sign -Wreorder -Wno-return-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings -Werror -diag-error:3846 > $<$<CXX_COMPILER_ID:MSVC>: /W4 >
# )


enable_testing()
include(CTest)

# add_executable(getrf.cpp.x getrf.cpp)
# add_test(NAME getrf.cpp.x COMMAND ./getrf.cpp.x)

# add_executable(geqrf.cpp.x geqrf.cpp)
# add_test(NAME geqrf.cpp.x COMMAND ./geqrf.cpp.x)

add_executable(potrf.cpp.x potrf.cpp)

if(ENABLE_CUDA)
    set_source_files_properties(potrf.cpp PROPERTIES LANGUAGE CUDA)
endif()

target_link_libraries(potrf.cpp.x PRIVATE multi-lapack Boost::boost)

add_test(NAME potrf.cpp.x COMMAND $<TARGET_FILE:potrf.cpp.x>)
