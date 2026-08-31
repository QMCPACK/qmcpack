cmake_minimum_required(VERSION 3.18)
project(
	boost-multi-adaptor-cuda-thrust-test
	VERSION 0.1
	LANGUAGES CXX CUDA
)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
	set(CMAKE_CUDA_ARCHITECTURES native)
endif()

if(MSVC)
	set(CMAKE_CUDA_RUNTIME_LIBRARY "None")  # "Shared"
endif()

find_package(CUDAToolkit REQUIRED)

# set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda --Werror=cross-execution-space-call -Xcudafe \"--display_error_number\"")

find_package(Boost CONFIG)  # tests require Boost.Core LightweightTest

include(CTest)
enable_testing()

# find_program(MEMORYCHECK_COMMAND valgrind)
# set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --error-exitcode=1")

set(TEST_SRCS
	array.cu
	broadcast_softmax.cu
	by_key.cu
	device_function.cu
	diffusion3d_kernel.cu
	fill_async.cu
	for_each.cu
	memory_resource.cu
	neighbor_list.cu
	restriction.cu
	run_reduce.cu
	set_identity_kernel.cu
	speed.cu
	speed_algo.cu
	transform.cu
	universal.cu
	vector.cu
)

foreach(TEST_FILE ${TEST_SRCS})
	set(TEST_EXE "${TEST_FILE}.x")
	add_executable(${TEST_EXE} ${TEST_FILE})
	if(ENABLE_CUDA OR DEFINED CXXCUDA)
		set_source_files_properties(${TEST_FILE} PROPERTIES LANGUAGE CUDA)
		# set_property(TARGET ${TEST_EXE} PROPERTY CUDA_STANDARD 17)
	endif()

	target_include_directories(${TEST_EXE}        PUBLIC ../../../../../include)

	target_compile_definitions(${TEST_EXE} PRIVATE "BOOST_PP_VARIADICS")
	target_compile_definitions(${TEST_EXE} PRIVATE ${Boost_DEFINITIONS})
	target_include_directories(${TEST_EXE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS})
	target_link_libraries(${TEST_EXE} PRIVATE ${Boost_LIBRARIES})
	target_link_libraries(${TEST_EXE} PRIVATE multi)
	target_link_libraries(${TEST_EXE} PRIVATE CUDA::cudart)
	target_link_directories(${TEST_EXE} PRIVATE ${Boost_LIBRARY_DIRS})
    if(ENABLE_CUDA)
	    add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
    endif()

	target_compile_options(${TEST_EXE} PRIVATE
		$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
			--diag-suppress=1394  # field of class type without a DLL interface used in a class with a DLL interface
			--extended-lambda
		#>
		$<$<CXX_COMPILER_ID:MSVC>:
		#	/GS      # Enables buffer security check (can help detect stack-based buffer overflows)
			$<$<NOT:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:
				/W4    # /volatile:iso /WX /Zc:wchar_t /Zc:forScope /Zc:inline
				/we4244  # treat 'conversion, possible loss of data' (narrowing) as an error
			>
			$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:
				--diag-suppress=1394  # field of class type without a DLL interface used in a class with a DLL interface
				-Xcompiler=/W4
				-Xcompiler=/we4100  # treat 'unreferenced formal parameter' as an error
				-Xcompiler=/we4244  # treat 'conversion, possible loss of data' (narrowing) as an error
				-Xcompiler=/wd4061  # enumerator 'cudaErrorSystemNotReady' in switch of enum 'cudaError' is not explicitly handled by a case label
				-Xcompiler=/wd4324  # structure was padded due to alignment specifier
				-Xcompiler=/wd4464
				-Xcompiler=/wd4664
				-Xcompiler=/wd4668
				-Xcompiler=/wd4296  # '<': expression is always false (from CUDA Toolkit's own cuda/std/__type_traits/is_signed.h)
			>
		>
	>
	)
endforeach()

set_tests_properties(speed.cu.x PROPERTIES RUN_SERIAL TRUE)

# target_compile_options(array.cu.x PRIVATE       $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:$<$<CXX_COMPILER_ID:MSVC>:-Xcompiler=/wd4668 -Xcompiler=/wd4664 -Xcompiler=/wd4464 > --extended-lambda >)
# target_compile_options(by_key.cu.x PRIVATE      $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:$<$<CXX_COMPILER_ID:MSVC>:-Xcompiler=/wd4668 -Xcompiler=/wd4664 -Xcompiler=/wd4464 > >)
# target_compile_options(run_reduce.cu.x PRIVATE  $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--extended-lambda >)
# target_compile_options(speed.cu.x PRIVATE       $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--extended-lambda >)

# target_compile_features(neighbor_list.cu.x PRIVATE cuda_std_20 cxx_std_20)
