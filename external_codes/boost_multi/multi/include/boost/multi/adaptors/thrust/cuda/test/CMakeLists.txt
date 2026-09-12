cmake_minimum_required(VERSION 3.16)

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
	find_package(CUDAToolkit REQUIRED)  # requires cmake 3.17
	# Thrust ships with the CUDA toolkit, but its CMake config lands in different places: apt puts it in
	# CMake's default search path (/usr/share/cmake/thrust), while NVIDIA toolkits/containers put it under
	# the toolkit (/usr/local/cuda/lib64/cmake/thrust) which is NOT searched by default. Hint at the toolkit,
	# and fall back to the (header-only) toolkit headers when no Thrust CMake package is installed at all.
	find_package(Thrust CONFIG QUIET HINTS "${CUDAToolkit_LIBRARY_DIR}/cmake")
	if(Thrust_FOUND)
		thrust_create_target(MyThrust)
	else()
		add_library(MyThrust INTERFACE)
		target_link_libraries(MyThrust INTERFACE CUDA::cudart)
		target_include_directories(MyThrust SYSTEM INTERFACE ${CUDAToolkit_INCLUDE_DIRS})
	endif()
endif()


find_package(Boost CONFIG REQUIRED)  # tests require Boost.Core LightweightTest

enable_testing()
include(CTest)

include_directories(${CMAKE_BINARY_DIR})

set(TEST_SRCS)
#	neighbor_list.cu)

foreach(TEST_FILE ${TEST_SRCS})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set(TEST_EXE "${TEST_FILE}.x")
		add_executable(${TEST_EXE} ${TEST_FILE})

		target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})

		# target_include_directories(${TEST_EXE}        PRIVATE ${PROJECT_SOURCE_DIR}/include)
		target_include_directories(${TEST_EXE}        PRIVATE ${CUDA_INCLUDE_DIRS})
	
		#target_compile_options(${TEST_EXE} PRIVATE --expt-relaxed-constexpr)

		target_include_directories(${TEST_EXE} SYSTEM PUBLIC ${_THRUST_INCLUDE_DIR})
		target_link_libraries(${TEST_EXE} PRIVATE MyThrust)

		target_link_libraries(${TEST_EXE} PRIVATE ${CUDA_LIBRARIES})
		target_link_libraries(${TEST_EXE} PRIVATE multi)
		target_link_libraries(${TEST_EXE} PRIVATE MyThrust)  # inherits Thrust usage requirements, incl. -Xcompiler=/wd4668 on MSVC+nvcc

		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		
		add_test(NAME ${TEST_EXE} COMMAND ./${TEST_EXE})
	endif()
endforeach()

# target_compile_options(neighbor_list.cu.x PRIVATE $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--extended-lambda>)
