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

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NCCL
  REQUIRED_VARS NCCL_LIB NCCL_INCLUDE_DIR
)

if(NCCL_FOUND AND NOT TARGET NCCL::NCCL)
  add_library(NCCL::NCCL INTERFACE IMPORTED)
  target_include_directories(NCCL::NCCL INTERFACE "${NCCL_INCLUDE_DIR}")
  target_link_libraries(NCCL::NCCL INTERFACE "${NCCL_LIB}")
endif()
