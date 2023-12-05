# -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
# Copyright 2018-2022 Alfredo A. Correa

cmake_minimum_required(VERSION 3.13)
# for reference: Ubuntu 20.04 ships cmake 3.16.3, Ubuntu 21.10 -> 3.18.4, Fedora 31 -> 3.17.4, (Gitlab CI) debian:stable-backports ships cmake 3.13.4, debian:testing -> 3.18.4

# to install this project: cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$HOME cmake --build . --config Release --target test --target install -- -j $(nproc)
# to use this project do: project("Your project") find_package(boost-multi CONFIG REQUIRED) add_executable(${PROJECT_NAME} src/your_main.cpp) target_link_libraries(${PROJECT_NAME} boost-multi::boost-multi)

message("CMake version: ${CMAKE_VERSION}")
message("current build directory: ${CMAKE_CURRENT_BINARY_DIR}")

project(
	multi
	VERSION 0.79.0
	DESCRIPTION "A header only C++ library that provides multidimensional array access to contiguous or regularly contiguous memory (or ranges)."
	HOMEPAGE_URL "https://gitlab.com/correaa/boost-multi"
	LANGUAGES CXX
)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
# 		$<$<AND:$<CXX_COMPILER_ID:GNU>,$<CUDA_COMPILER_ID:NVIDIA>>:
# https://youtrack.jetbrains.com/issue/CPP-25608
target_compile_features(${PROJECT_NAME} INTERFACE $<$<NOT:$<CUDA_COMPILER_ID:NVIDIA>>:cxx_std_17>)
target_compile_options (${PROJECT_NAME} INTERFACE       $<$<CUDA_COMPILER_ID:NVIDIA>: -std=c++17>)

target_compile_options (
	${PROJECT_NAME}
	INTERFACE
		$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
			-Xcompiler=-Werror,-Wall,-Wextra,-Wcast-align,-Wcast-qual,-Wno-double-promotion,-Wduplicated-branches,-Wduplicated-cond,-Wformat-truncation,-Wformat=2,-Wlogical-op,-Wmisleading-indentation,-Wno-missing-include-dirs,-Wnon-virtual-dtor,-Wno-missing-declarations,-Wnon-virtual-dtor,-Wnull-dereference,-Woverloaded-virtual,-Wpointer-arith,-Wno-redundant-decls,-Wno-shadow,-Wno-switch-enum,-Wno-unknown-pragmas,-Wtrampolines,-Wuninitialized,-Wunused,-Wunused-but-set-variable,-Wunused-result,-Wno-zero-as-null-pointer-constant
			--expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe=--display_error_number -Xcudafe=--diag_error=incompatible_assignment_operands -Xcudafe=--diag_error=returning_ptr_to_local_variable -Xcudafe=--diag_error=subscript_out_of_range -Xcudafe=--diag_error=used_before_set -Xcudafe=--diag_error=undefined_preproc_id -Xcudafe=--diag_error=implicit_func_decl -Xcudafe=--diag_error=implicit_return_from_non_void_function -Xcudafe=--diag_error=missing_type_specifier
		>
		$<$<OR:$<CXX_COMPILER_ID:PGI>,$<CXX_COMPILER_ID:NVHPC>>:  # EDG diagnostics list: http://www.ssl.berkeley.edu/~jimm/grizzly_docs/SSL/opt/intel/cc/9.0/lib/locale/en_US/mcpcom.msg
        	--display_error_number --diag_error=incompatible_assignment_operands,returning_ptr_to_local_variable,subscript_out_of_range,used_before_set,undefined_preproc_id,implicit_func_decl,implicit_return_from_non_void_function,missing_type_specifier
		>
)

install(
	TARGETS ${PROJECT_NAME}
	EXPORT ${PROJECT_NAME}_Targets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${PROJECT_NAME}ConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

install(
	EXPORT ${PROJECT_NAME}_Targets
	FILE ${PROJECT_NAME}Targets.cmake
	NAMESPACE boost::multi::
	DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

if (NOT MSVC)
	install(FILES     ${PROJECT_SOURCE_DIR}/include/multi/array_ref.hpp ${PROJECT_SOURCE_DIR}/include/multi/array.hpp DESTINATION include/multi)
	install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/multi/detail                                                      DESTINATION include/multi)
	install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/multi/adaptors                                                    DESTINATION include/multi)
endif()

#find_program(MEMORYCHECK_COMMAND valgrind)
#set(VALGRIND_COMMAND_OPTIONS    "-q --tool=memcheck --leak-check=yes --num-callers=51 --trace-children=yes --leak-check=full --track-origins=yes --gen-suppressions=all")
#set(MEMORYCHECK_COMMAND_OPTIONS "-q --tool=memcheck --leak-check=yes --num-callers=52 --trace-children=yes --leak-check=full --track-origins=yes --gen-suppressions=all")  # must go before `include(CTest)`
#set(MEMORYCHECK_SUPPRESSIONS_FILE "${PROJECT_SOURCE_DIR}/.valgrind-suppressions")  # must go before `include(CTest)`

include(CTest)

find_package(Boost COMPONENTS unit_test_framework)

if(Boost_FOUND)
	enable_testing()
	add_subdirectory(test)

	add_subdirectory(include/multi/adaptors/blas  )
	add_subdirectory(include/multi/adaptors/fftw  )
	add_subdirectory(include/multi/adaptors/cuda  )
	add_subdirectory(include/multi/adaptors/thrust)
else()
    message (WARNING "Cannot find Boost, library will not be tested. If you want this feature install Boost.Test, for example please run:\n    sudo apt install libboost-test-dev")
endif()

