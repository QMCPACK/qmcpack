find_path(NCCL_INCLUDE_DIR NAMES nccl.h
  HINTS
    "${NCCL_ROOT}/include"
    "$ENV{NCCL_ROOT}/include"
    "${CUDAToolkit_INCLUDE_DIRS}"
    "${CUDA_TOOLKIT_ROOT_DIR}/include"
    "$ENV{CUDA_PATH}/include"
)

find_library(NCCL_LIB nccl
  HINTS
    "${NCCL_ROOT}"
    "$ENV{NCCL_ROOT}"
    "${CUDAToolkit_LIBRARY_DIR}"
    "${CUDA_TOOLKIT_ROOT_DIR}"
    "$ENV{CUDA_PATH}"
  PATH_SUFFIXES lib64 lib
)

if(NCCL_INCLUDE_DIR AND EXISTS "${NCCL_INCLUDE_DIR}/nccl.h")
  file(STRINGS "${NCCL_INCLUDE_DIR}/nccl.h" _nccl_version_lines
       REGEX "^[ \t]*#[ \t]*define[ \t]+NCCL_(MAJOR|MINOR|PATCH)[ \t]+[0-9]+")
  foreach(_line IN LISTS _nccl_version_lines)
    if(_line MATCHES "NCCL_MAJOR[ \t]+([0-9]+)")
      set(NCCL_VERSION_MAJOR "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "NCCL_MINOR[ \t]+([0-9]+)")
      set(NCCL_VERSION_MINOR "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "NCCL_PATCH[ \t]+([0-9]+)")
      set(NCCL_VERSION_PATCH "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(DEFINED NCCL_VERSION_MAJOR AND DEFINED NCCL_VERSION_MINOR AND DEFINED NCCL_VERSION_PATCH)
    set(NCCL_VERSION "${NCCL_VERSION_MAJOR}.${NCCL_VERSION_MINOR}.${NCCL_VERSION_PATCH}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NCCL
  REQUIRED_VARS NCCL_LIB NCCL_INCLUDE_DIR
  VERSION_VAR NCCL_VERSION
)

if(NCCL_FOUND AND NOT TARGET NCCL::NCCL)
  add_library(NCCL::NCCL INTERFACE IMPORTED)
  target_include_directories(NCCL::NCCL INTERFACE "${NCCL_INCLUDE_DIR}")
  target_link_libraries(NCCL::NCCL INTERFACE "${NCCL_LIB}")
endif()
