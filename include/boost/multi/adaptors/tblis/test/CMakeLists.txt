cmake_minimum_required(VERSION 3.16)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(
	boost-multi-adaptors-blis-test
	VERSION 0.1
	LANGUAGES CXX
)

option(ENABLE_TBLIS "Enable TBLIS tests" OFF)

if(NOT ENABLE_TBLIS)
	message(WARNING "TBLIS disabled, will not be tested")
	return()
else()
	message(WARNING "TBLIS enabled, will be installed and tested")

	include(FetchContent)

	enable_testing()
	include(CTest)

	# include_directories(${CMAKE_BINARY_DIR})

	#include_directories(/home/correaa/include /home/correaa/tblis/src/external/tci)
	#link_libraries(/home/correaa/lib/libtblis.so)

	# set(OpenMP_C_FLAGS "-Xclang;-fopenmp" CACHE STRING "" FORCE)
	# set(OpenMP_CXX_FLAGS "-Xclang;-fopenmp" CACHE STRING "" FORCE)

	# set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "" FORCE)
	# set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)

	include(ExternalProject)

	# Where all externals are installed
	set(EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external)

	# --- CRITICAL: make directories exist at configure time ---
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_DIR})
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_DIR}/include)
	file(MAKE_DIRECTORY ${EXTERNAL_INSTALL_DIR}/lib)

	# ------------------------------------------------------------------------------
	# TBLIS external project
	# ------------------------------------------------------------------------------

	ExternalProject_Add(tblis_ep
		# URL https://github.com/MatthewsResearchGroup/tblis/archive/refs/tags/v2.0-beta2.tar.gz
		# PATCH_COMMAND git submodule update --init --recursive
		GIT_REPOSITORY https://github.com/devinamatthews/tblis.git
		GIT_TAG master  # v2.0-beta2  # pin to a commit/tag in real projects
		# GIT_SHALLOW TRUE

		PREFIX ${CMAKE_BINARY_DIR}/tblis

		CMAKE_ARGS
			-DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_DIR}
			-DCMAKE_BUILD_TYPE=Release
			-DCMAKE_C_COMPILER=gcc

		BUILD_BYPRODUCTS
			${EXTERNAL_INSTALL_DIR}/lib/libtblis${CMAKE_SHARED_LIBRARY_SUFFIX}

		BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${PROCESSOR_COUNT}
	)

	# Paths *after* installation
	set(TBLIS_INCLUDE_DIR ${EXTERNAL_INSTALL_DIR}/include)
	set(TBLIS_LIBRARY
		${EXTERNAL_INSTALL_DIR}/lib/libtblis${CMAKE_SHARED_LIBRARY_SUFFIX}
	)

	add_library(tblis::tblis SHARED IMPORTED GLOBAL)
	#target_compile_features(tblis::tblis PUBLIC cxx_std_20)

	set_target_properties(tblis::tblis PROPERTIES
		IMPORTED_LOCATION ${TBLIS_LIBRARY}
		INTERFACE_INCLUDE_DIRECTORIES
			$<BUILD_INTERFACE:${TBLIS_INCLUDE_DIR}>
	)

	# Ensure build order
	add_dependencies(tblis::tblis tblis_ep)

	#target_link_libraries(tblis::tblis PRIVATE $<$<CXX_COMPILER_ID:GNU>:atomic>)  # for macos

	# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
	set(TEST_SRCS tensor.cpp)

	foreach(TEST_FILE ${TEST_SRCS})
		set(TEST_EXE "${TEST_FILE}.x")
		add_executable(${TEST_EXE} ${TEST_FILE})
		if(ENABLE_CUDA OR DEFINED CXXCUDA)
			set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
			# target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
		endif()
		# target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17)
		# target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
		target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
		target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})

		# target_link_libraries(${TEST_EXE} PRIVATE multi)
		target_link_libraries(${TEST_EXE} PRIVATE tblis::tblis multi)

		target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
		target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})

		if(NOT ENABLE_CUDA)
			target_compile_options(
				${TEST_EXE}
				PRIVATE -Werror
						-fno-common
						$<$<CXX_COMPILER_ID:GNU>:
						-Wformat-truncation
						-fstack-usage> # -Wconversion
						$<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
						-Wmove>
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
						#/W4
						>
			)
		endif()
		add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
	endforeach()
endif()

