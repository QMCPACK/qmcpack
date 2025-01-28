#--------------------------------------------------------------------
# QMC_GPU option
#--------------------------------------------------------------------

# Guess QMC_GPU_DEFAULT if QMC_GPU_ARCHS is provided. Only use the first entry if there are more.
if(QMC_GPU_ARCHS)
  list(GET QMC_GPU_ARCHS 0 GPU_ARCH_0)
  if(GPU_ARCH_0 MATCHES "^sm_")
    set(QMC_GPU_DEFAULT "openmp;cuda")
  elseif(GPU_ARCH_0 MATCHES "^gfx")
    set(QMC_GPU_DEFAULT "openmp;hip")
  elseif(GPU_ARCH_0 MATCHES "^intel_gpu_")
    set(QMC_GPU_DEFAULT "openmp;sycl")
  else()
    set(QMC_GPU_DEFAULT "")
  endif()
endif()

set(VALID_QMC_GPU_FEATURES "openmp" "cuda" "hip" "sycl")
set(QMC_GPU ${QMC_GPU_DEFAULT} CACHE STRING "Semicolon-separated list of GPU features to enable (openmp,cuda,hip,sycl)")
set_property(CACHE QMC_GPU PROPERTY STRINGS ${VALID_QMC_GPU_FEATURES})

set(ENABLE_OFFLOAD OFF)
set(ENABLE_CUDA OFF)
set(QMC_CUDA2HIP OFF)
set(ENABLE_SYCL OFF)

# Perform QMC_GPU option check
foreach(GPU_FEATURE IN LISTS QMC_GPU)
  # verify the entry
  if(NOT GPU_FEATURE IN_LIST VALID_QMC_GPU_FEATURES)
    message(FATAL_ERROR "Invalid QMC_GPU selection \"${GPU_FEATURE}\", valid values are \"${VALID_QMC_GPU_FEATURES}\"")
  endif()

  if(GPU_FEATURE STREQUAL "openmp")
    set(ENABLE_OFFLOAD ON)
  endif()
  if(GPU_FEATURE STREQUAL "cuda")
    set(ENABLE_CUDA ON)
  endif()
  if(GPU_FEATURE STREQUAL "hip")
    set(ENABLE_CUDA ON)
    set(QMC_CUDA2HIP ON)
  endif()
  if(GPU_FEATURE STREQUAL "sycl")
    set(ENABLE_SYCL ON)
  endif()
endforeach()

if("cuda" IN_LIST QMC_GPU AND "hip" IN_LIST QMC_GPU OR
   "cuda" IN_LIST QMC_GPU AND "sycl" IN_LIST QMC_GPU OR
   "hip" IN_LIST QMC_GPU AND "sycl" IN_LIST QMC_GPU)
    message(FATAL_ERROR "Invalid QMC_GPU selection\"${QMC_GPU}\", \"cuda\", \"hip\" and \"sycl\" don't work together. Select one!")
endif()

if(QMC_GPU)
  message(STATUS "Enable GPU features QMC_GPU=${QMC_GPU}")
endif()

if(QMC_CUDA2HIP)
  set(ENABLE_ROCM ON)
else()
  set(ENABLE_ROCM OFF)
endif()
