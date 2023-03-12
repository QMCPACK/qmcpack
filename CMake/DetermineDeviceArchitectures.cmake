# derive AMD GPU architectures from CMAKE_HIP_ARCHITECTURES or detect them using amdgpu-arch (LLVM15+)
function(detectAMDGPU)
  if(CMAKE_HIP_ARCHITECTURES)
    set(QMC_GPU_ARCHS_DETECTED ${CMAKE_HIP_ARCHITECTURES})
  else()
    find_program(AMDGPU_ARCH_EXE amdgpu-arch)
    if(AMDGPU_ARCH_EXE)
      execute_process(COMMAND ${AMDGPU_ARCH_EXE} OUTPUT_VARIABLE AMD_GPU_ARCH)
      string(REGEX REPLACE "\n" ";" AMD_GPU_ARCH "${AMD_GPU_ARCH}")
      if(AMD_GPU_ARCH)
        list(APPEND QMC_GPU_ARCHS_DETECTED ${AMD_GPU_ARCH})
      else()
        message("amdgpu-arch didn't find AMD GPUs. Cannot auto-detect AMD GPU architectures.")
      endif()
    else()
      message("amdgpu-arch not avaible. Cannot auto-detect AMD GPU architectures.")
    endif()
  endif()
  set(QMC_GPU_ARCHS
      ${QMC_GPU_ARCHS_DETECTED}
      PARENT_SCOPE)
endfunction()

# check QMC_GPU_ARCHS and CMAKE_HIP_ARCHITECTURES consistency
function(verifyAMDGPUconsistency)
  if(CMAKE_HIP_ARCHITECTURES AND NOT CMAKE_HIP_ARCHITECTURES STREQUAL QMC_GPU_ARCHS)
    message(FATAL_ERROR "CMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES} doesn't match ${QMC_GPU_ARCHS}.")
  endif()
endfunction()

# derive NVIDIA GPU architectures from CMAKE_CUDA_ARCHITECTURES or detect them using nvptx-arch (LLVM16+)
function(detectNVIDIAGPU)
  if(CMAKE_CUDA_ARCHITECTURES)
    set(QMC_GPU_ARCHS_DETECTED)
    foreach(CUDA_ARCH_NUM IN LISTS CMAKE_CUDA_ARCHITECTURES)
      list(APPEND QMC_GPU_ARCHS_DETECTED "sm_${CUDA_ARCH_NUM}")
    endforeach()
  else()
    find_program(NVPTX_ARCH_EXE nvptx-arch)
    if(NVPTX_ARCH_EXE)
      execute_process(COMMAND ${NVPTX_ARCH_EXE} OUTPUT_VARIABLE NVIDIA_GPU_ARCH)
      string(REGEX REPLACE "\n" ";" NVIDIA_GPU_ARCH "${NVIDIA_GPU_ARCH}")
      if(NVIDIA_GPU_ARCH)
        list(APPEND QMC_GPU_ARCHS_DETECTED ${NVIDIA_GPU_ARCH})
      else()
        message("nvptx-arch didn't find NVIDIA GPUs. Cannot auto-detect NVIDIA GPU architectures.")
      endif()
    else()
      message("nvptx-arch not avaible. Cannot auto-detect NVIDIA GPU architectures.")
    endif()
  endif()
  set(QMC_GPU_ARCHS
      ${QMC_GPU_ARCHS_DETECTED}
      PARENT_SCOPE)
endfunction()

# check QMC_GPU_ARCHS and CMAKE_CUDA_ARCHITECTURES consistency
function(verifyNVIDIAGPUconsistency)
  string(REPLACE "sm_" "" CUDA_ARCH_NUMBERS "${QMC_GPU_ARCHS}")
  if(CMAKE_CUDA_ARCHITECTURES AND NOT CMAKE_CUDA_ARCHITECTURES STREQUAL CUDA_ARCH_NUMBERS)
    message(FATAL_ERROR "CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES} doesn't match ${QMC_GPU_ARCHS}.")
  endif()
endfunction()

# auto detect QMC_GPU_ARCHS if not set by user and GPU features are enabled.
if(NOT QMC_GPU_ARCHS AND ENABLE_CUDA)
  if(QMC_CUDA2HIP)
    detectAMDGPU()
  else()
    detectNVIDIAGPU()
  endif()

  if(NOT QMC_GPU_ARCHS)
    message(
      FATAL_ERROR
        "QMC_GPU_ARCHS was neither set nor derivable and auto-detection failed. "
        "QMC_GPU_ARCHS is required to be set for enabling CUDA/HIP features. "
        "For example, set sm_80 for NVIDIA A100 GPU and gfx90a for AMD MI200 series GPUs.")
  endif()
endif()

list(REMOVE_DUPLICATES QMC_GPU_ARCHS)

# make sure QMC_GPU_ARCHS is consistent with CMAKE_HIP_ARCHITECTURES or CMAKE_CUDA_ARCHITECTURES.
if(ENABLE_CUDA)
  if(QMC_CUDA2HIP)
    verifyAMDGPUconsistency()
  else()
    verifyNVIDIAGPUconsistency()
  endif()
endif()

set(QMC_GPU_ARCHS
    ${QMC_GPU_ARCHS}
    CACHE STRING "Accelerator device architectures" FORCE)
