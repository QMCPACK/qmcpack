cmake_minimum_required(VERSION 3.22)

project(
	boost-multi-adaptors-vkfft-test
	VERSION 0.1
	LANGUAGES CXX
)

# multi::vkfft (VkFFT headers + backend libs + VKFFT_BACKEND) comes from the parent
# CMakeLists.  Pick the source language from the enabled GPU toolchain.
if(ENABLE_CUDA OR DEFINED CXXCUDA)
	set(_vkfft_lang CUDA)
	if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES native)
	endif()
elseif(ENABLE_HIP)
	set(_vkfft_lang HIP)
else()
	return()
endif()
enable_language(${_vkfft_lang})

find_package(Boost CONFIG REQUIRED)                    # <boost/core/lightweight_test.hpp>

# FFTW is the CPU cross-check reference for the tests.
find_package(PkgConfig REQUIRED)  # install with `winget install pkgconf.pkgconf`
pkg_search_module(FFTW REQUIRED IMPORTED_TARGET fftw3)  # sudo apt install libfftw3-dev

add_executable(vkfft.cpp.x vkfft.cpp)
set_source_files_properties(vkfft.cpp PROPERTIES LANGUAGE ${_vkfft_lang})
target_link_libraries(vkfft.cpp.x PRIVATE multi::vkfft Boost::headers PkgConfig::FFTW)
add_test(NAME vkfft.cpp.x COMMAND vkfft.cpp.x)
