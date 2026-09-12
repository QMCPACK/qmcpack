cmake_minimum_required(VERSION 3.22)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(
	boost-multi-adaptors-vkfft
	VERSION 0.1
	LANGUAGES CXX
)

# multi::vkfft - the boost::multi VkFFT adaptor plus everything a translation unit
# needs to use it: the VkFFT headers, the backend's runtime-compile libraries, and
# the VKFFT_BACKEND selection.  (Not yet part of the installed `multi` export, same
# as multi_fftw3 / multi-lapack / multi-blas.)
add_library(multi_vkfft INTERFACE)
add_library(multi::vkfft ALIAS multi_vkfft)
target_link_libraries(multi_vkfft INTERFACE multi)

# VkFFT headers: prefer an installed copy (FindVkFFT.cmake), otherwise fetch them.
find_package(VkFFT QUIET)
if(NOT VkFFT_FOUND)
	include(FetchContent)
	FetchContent_Declare(
		vkfft
		GIT_REPOSITORY https://github.com/DTolm/VkFFT.git
		GIT_TAG        v1.3.4
		GIT_SHALLOW    TRUE
		SOURCE_SUBDIR  do-not-configure  # headers only; skip VkFFT's own CMakeLists
	)
	FetchContent_MakeAvailable(vkfft)
	add_library(VkFFT::VkFFT INTERFACE IMPORTED)
	set_target_properties(VkFFT::VkFFT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${vkfft_SOURCE_DIR}/vkFFT")
	message(STATUS "boost-multi vkfft: using fetched VkFFT (v1.3.4); install VkFFT and set -DVkFFT_ROOT to skip the download")
endif()
target_link_libraries(multi_vkfft INTERFACE VkFFT::VkFFT)

# set_property(TARGET multi_vkfft APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
# 	$<TARGET_PROPERTY:VkFFT::VkFFT,INTERFACE_INCLUDE_DIRECTORIES>)

# Backend selection.  VKFFT_BACKEND is a compile-time choice; set it (and pull the
# runtime compiler) from whichever GPU toolchain is enabled.  The header also
# auto-detects from __CUDACC__ / __HIPCC__, so a plain nvcc/hipcc build still works
# without this define.
if(ENABLE_CUDA OR DEFINED CXXCUDA)
	find_package(CUDAToolkit REQUIRED)
	target_link_libraries(multi_vkfft INTERFACE CUDA::cudart CUDA::cuda_driver CUDA::nvrtc)
	target_compile_definitions(multi_vkfft INTERFACE VKFFT_BACKEND=1)
elseif(ENABLE_HIP)
	enable_language(HIP)
	find_package(hip REQUIRED)
	target_link_libraries(multi_vkfft INTERFACE hip::host)

	# VkFFT's HIP backend generates kernels at runtime with HIPRTC.  On ROCm < 7 the
	# HIPRTC symbols live inside libamdhip64 (already pulled in by hip::host); on
	# ROCm >= 7 HIPRTC is a separate library with its own package.
	find_package(hiprtc QUIET)
	if(TARGET hiprtc::hiprtc)
		target_link_libraries(multi_vkfft INTERFACE hiprtc::hiprtc)
	else()
		find_library(VkFFT_HIPRTC_LIBRARY NAMES hiprtc HINTS ${ROCM_PATH} ENV ROCM_PATH /opt/rocm PATH_SUFFIXES lib lib64)
		if(VkFFT_HIPRTC_LIBRARY)
			target_link_libraries(multi_vkfft INTERFACE ${VkFFT_HIPRTC_LIBRARY})
		endif()
		mark_as_advanced(VkFFT_HIPRTC_LIBRARY)
	endif()

	target_compile_definitions(multi_vkfft INTERFACE VKFFT_BACKEND=2)
endif()

add_subdirectory(test)  # the test CMakeLists returns early when no GPU toolchain is enabled
