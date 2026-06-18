cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(
	boost-multi-adaptors-blas-test
	VERSION 0.1
	LANGUAGES CXX
)

# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT
   CMAKE_CXX_COMPILER_ID
   STREQUAL
   "PGI"
)
	find_package(Boost CONFIG)
else()
	find_package(Boost CONFIG)
	# find_package(Boost REQUIRED) # cmake cannot detect this component with pgi compiler
endif()

if(ENABLE_CUDA OR DEFINED CXXCUDA)
	enable_language(CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()

	# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr --extended-lambda")
	# find_package(CUDA QUIET)
	find_package(CUDAToolkit REQUIRED COMPONENTS cuFFT)
	# if(CUDA_FOUND)
	#   message("CUDA found")
	#   include_directories(${CUDA_INCLUDE_DIRS})
	#   # set(INQ_CUDA_LIBRARIES ${CUDA_CUFFT_LIBRARIES} ${CUDA_cusolver_LIBRARY} ${CUDA_CUBLAS_LIBRARIES})
	#   # link_libraries(${INQ_CUDA_LIBRARIES})
	# else()
	#   message("CUDA not found")
	# endif()
endif()

# FFTW has to go before blas to avoid unscrupulous (i.e. MKL) blas implementations that include FFTW and don't implement it properly
find_package(PkgConfig)

if(NOT PKG_CONFIG_FOUND)
	message(WARNING "pkgconfig  not found (for FFTW). cuFFT tests needs FFTW.")
	return()
endif()

pkg_search_module(  # install fftw if necessary: sudo apt install libfftw3-dev 
	FFTW
	REQUIRED
	fftw3
	IMPORTED_TARGET
)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

enable_testing()
list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure") # needs cmake 3.17
include(CTest)

include_directories(${CMAKE_BINARY_DIR})

# file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
set(TEST_SRCS cufft.cpp)

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
	target_link_libraries(${TEST_EXE} PRIVATE CUDA::cufft)

	add_test(NAME ${TEST_EXE} COMMAND $<TARGET_FILE:${TEST_EXE}>)
endforeach()
