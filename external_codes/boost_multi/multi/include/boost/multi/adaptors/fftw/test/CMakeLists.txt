cmake_minimum_required(VERSION 3.10)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_VERBOSE_MAKEFILE ON)

#project(
#   boost-multi-adaptors-fftw-test
#   VERSION 0.1
#   LANGUAGES CXX
#)

if((NOT
	CMAKE_CXX_COMPILER_ID
	STREQUAL
	"PGI"
   )
   AND (NOT
		DART_COMPILER_NAME
		STREQUAL
		"nvcc"
	   )
)
#   find_package(Boost REQUIRED COMPONENTS unit_test_framework timer)
else()
	link_libraries("-lboost_unit_test_framework") # -lboost_timer")
endif()

# find_package(PkgConfig REQUIRED)
# pkg_search_module(
#   FFTW
#   REQUIRED
#   fftw3
#   IMPORTED_TARGET
# )

# include_directories(PkgConfig::FFTW)
# link_libraries(PkgConfig::FFTW)
# cmake_print_properties(TARGETS PkgConfig::FFTW PROPERTIES INTERFACE_LINK_LIBRARIES INTERFACE_INCLUDE_DIRECTORIES)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
endif()

#find_package(CUDA QUIET)

#if(CUDA_FOUND)
#   message("CUDA found")
#   include_directories(${CUDA_INCLUDE_DIRS})
#else()
#   message("CUDA not found")
#endif()

enable_testing()

include(CTest)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS
  combinations.cpp
# copy.cpp
  core.cpp
  moved.cpp
  shift.cpp
  so_shift.cpp
# thrust.cpp
  transpose.cpp
  transpose_square.cpp
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
	#   set_property(TARGET ${TEST_EXE} PROPERTY "${CUDA_ARCH_LIST}")
		target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
	endif()

	target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
	if((NOT ENABLE_CUDA)
	   AND (NOT
			DART_COMPILER_NAME
			STREQUAL
			"nvcc"
		   )
	)
        if(NOT DEFINED ENABLE_CIRCLE)
		  target_compile_options(
			  ${TEST_EXE}
			  PRIVATE -Werror
					  -Wall
					  -Wextra
					  $<$<CXX_COMPILER_ID:GNU>:
					  -fno-common
					  -Wpedantic
					  -Wformat-truncation
					  -Wparentheses
					  -fstack-usage> # -Wconversion
					  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
					  -fno-common
					  -Wpedantic
					  -Wmove
					  -Wno-error=\#warnings
					  >
					  $<$<CXX_COMPILER_ID:Intel>:
					  -wd161
					  -diag-disable=remark
					  -Warray-bounds
					  -Wchar-subscripts
					  -Wcomment
					  -Wenum-compare
					  -Wformat
					  -Wuninitialized
					  -Wmaybe-uninitialized
					  -Wmain
					  -Wnarrowing
					  -Wnonnull
					  -Wparentheses
					  -Wpointer-sign
					  -Wreorder
					  -Wno-return-type
					  -Wsign-compare
					  -Wsequence-point
					  -Wtrigraphs
					  -Wunused-function
					  -Wunused-but-set-variable
					  -Wunused-variable
					  -Wwrite-strings
					  -Werror
					  -diag-error:3846
					  >
					  $<$<CXX_COMPILER_ID:MSVC>:
					  /W4>
		  )
        else()
          target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_TEST_DYN_LINK=1")
        endif()
	endif()

	target_link_libraries(${TEST_EXE} PRIVATE multi PkgConfig::FFTW)

	add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
endforeach()
