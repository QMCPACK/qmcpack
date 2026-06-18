cmake_minimum_required(VERSION 3.17)

project(
	boost-multi-adaptors-cutensor-test
	VERSION 0.1
	LANGUAGES CXX
)

find_package(Boost CONFIG)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()

	# find_package(cuTENSOR REQUIRED)

	find_package(CUDAToolkit REQUIRED COMPONENTS cutensor)

	# add_library(cutensor::cutensor UNKNOWN IMPORTED)

	# set_target_properties(cutensor::cutensor PROPERTIES
	# 	IMPORTED_LOCATION "${CUDAToolkit_LIBRARY_DIR}/libcutensor.so"
	# 	INTERFACE_INCLUDE_DIRECTORIES "${CUDAToolkit_INCLUDE_DIRS}"
	# )

endif()

enable_testing()

include(CTest)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS cutensor.cu)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		# target_compile_options(${TEST_EXE} PRIVATE -std=c++17)
	endif()
	# target_compile_features   (${TEST_EXE} PUBLIC cxx_std_17)

	target_compile_options(
		${TEST_EXE} PRIVATE
		$<$<NOT:$<CXX_COMPILER_ID:Clang>>:-Wno-error=terminate>
	)

	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} PRIVATE ${Boost_INCLUDE_DIRS})

	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
	target_link_libraries(${TEST_EXE} PRIVATE multi)
	target_link_libraries(${TEST_EXE} PRIVATE cuTENSOR::cuTENSOR CUDA::cudart)

	add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
endforeach()
