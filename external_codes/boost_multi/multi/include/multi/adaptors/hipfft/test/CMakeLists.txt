cmake_minimum_required(VERSION 3.11)

set(CMAKE_VERBOSE_MAKEFILE ON)

project(
	boost-multi-adaptors-hipfft-test
	VERSION 0.1
	LANGUAGES CXX
)

enable_language(HIP)
find_package(hipFFT REQUIRED)

# FFTW has to go before blas to avoid unscrupulous (i.e. MKL) blas implementations that include FFTW and don't implement it properly
find_package(PkgConfig REQUIRED)
pkg_search_module(
	FFTW
	REQUIRED
	fftw3
	IMPORTED_TARGET
)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

enable_testing()
include(CTest)

# include_directories(${CMAKE_BINARY_DIR})

find_package(Boost COMPONENTS unit_test_framework)

add_executable(hipfft.cpp.x hipfft.cpp)
set_source_files_properties(hipfft.cpp PROPERTIES LANGUAGE HIP)

# target_link_libraries(${TEST_EXE} PRIVATE multi)
target_link_libraries(hipfft.cpp.x PRIVATE hip::hipfft Boost::unit_test_framework)

add_test(NAME hipfft.cpp.x COMMAND hipfft.cpp.x)
