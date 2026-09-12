# FindVkFFT.cmake - locate the header-only VkFFT library (github.com/DTolm/VkFFT)
#
# VkFFT ships no pkg-config file and no consumer-facing package config, so this
# module fills that gap, mirroring adaptors/cutensor/cmake/FindcuTENSOR.cmake.
#
# Result variables:
#   VkFFT_FOUND        - true if vkFFT.h was located
#   VkFFT_INCLUDE_DIR  - the directory containing vkFFT.h (and its vkFFT/ subtree)
#
# Imported target:
#   VkFFT::VkFFT       - INTERFACE target carrying VkFFT_INCLUDE_DIR
#
# Hint with -DVkFFT_ROOT=<prefix> (or the VkFFT_ROOT environment variable).

find_path(
	VkFFT_INCLUDE_DIR
	NAMES vkFFT.h
	HINTS ${VkFFT_ROOT} ENV VkFFT_ROOT
	PATH_SUFFIXES vkFFT VkFFT include include/vkFFT
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VkFFT REQUIRED_VARS VkFFT_INCLUDE_DIR)

if(VkFFT_FOUND AND NOT TARGET VkFFT::VkFFT)
	add_library(VkFFT::VkFFT INTERFACE IMPORTED)
	set_target_properties(VkFFT::VkFFT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${VkFFT_INCLUDE_DIR}")
endif()

mark_as_advanced(VkFFT_INCLUDE_DIR)
