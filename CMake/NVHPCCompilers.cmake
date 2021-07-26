# Enable OpenMP
# If just -mp is specified, OMP_NUM_THREADS must be set in order to run in parallel
# Specifying 'allcores' will run on all cores if OMP_NUM_THREADS is not set (which seems
#  to be the default for other OpenMP implementations)
if(QMC_OMP)
  set(ENABLE_OPENMP 1)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mp=allcores")
  if(ENABLE_OFFLOAD AND NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    message(WARNING "QMCPACK OpenMP offload is not ready for NVIDIA HPC compiler.")
    if(NOT DEFINED OFFLOAD_ARCH)
      message(FATAL_ERROR "NVIDIA HPC compiler requires -gpu=ccXX option set based on the target GPU architecture! "
                          "Please add -DOFFLOAD_ARCH=ccXX to cmake. For example, cc70 is for Volta.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mp=gpu")
    set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-gpu=${OFFLOAD_ARCH}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mp=allcores")
  endif()
endif(QMC_OMP)

add_definitions(-Drestrict=__restrict__)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__forceinline=inline")

# Suppress compile warnings
# 177 variable "XX" was declared but never referenced
# 550 variable "XX" was set but never used
# 612 overloaded virtual function "AA" is only partially overridden in class "BB"
# 998 function "AA" is hidden by "BB" -- virtual function override intended?
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --display_error_number --diag_suppress 177 --diag_suppress 550")
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} --display_error_number --diag_suppress 177 --diag_suppress 550 --diag_suppress 612 --diag_suppress 998"
)

# Set extra optimization specific flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast")

# Setting this to 'OFF' adds the -A flag, which enforces strict standard compliance
#  and causes the compilation to fail with some GNU header files
set(CMAKE_CXX_EXTENSIONS ON)

# Add static flags if necessary
if(QMC_BUILD_STATIC)
  set(CMAKE_CXX_LINK_FLAGS " -Bstatic")
endif(QMC_BUILD_STATIC)
