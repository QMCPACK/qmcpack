# Check compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
  message(FATAL_ERROR "Requires GCC 9.0 or higher ")
endif()

# Enable OpenMP
if(QMC_OMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

  if(ENABLE_OFFLOAD)
    message(WARNING "QMCPACK OpenMP offload is not ready for GCC compiler.")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
      message(WARNING "GCC OpenMP offload feature requires 12.0 or higher.")
    endif()

    if(QMC_CUDA2HIP)
      set(OFFLOAD_TARGET_DEFAULT "amdgcn-amdhsa")
    else()
      set(OFFLOAD_TARGET_DEFAULT "nvptx-none")
    endif()
    set(OFFLOAD_TARGET
        ${OFFLOAD_TARGET_DEFAULT}
        CACHE STRING "Offload target architecture")
    set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-foffload=${OFFLOAD_TARGET} -foffload-options=\"-lm -latomic\"")

    if(NOT DEFINED OFFLOAD_ARCH AND DEFINED QMC_GPU_ARCHS)
      list(LENGTH QMC_GPU_ARCHS QMC_GPU_ARCH_COUNT)
      if(QMC_GPU_ARCH_COUNT EQUAL "1")
        set(OFFLOAD_ARCH ${QMC_GPU_ARCHS})
      else()
        message(
          FATAL_ERROR
            "GCC does not yet support offload to multiple architectures! "
            "Deriving OFFLOAD_ARCH from QMC_GPU_ARCHS failed. "
            "Please keep only one entry in QMC_GPU_ARCHS or set OFFLOAD_ARCH.")
      endif()
    endif()

    if(DEFINED OFFLOAD_ARCH)
      if(OFFLOAD_TARGET MATCHES "amdgcn-amdhsa")
        set(OPENMP_OFFLOAD_COMPILE_OPTIONS
            "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -foffload-options=${OFFLOAD_TARGET}=\"-march=${OFFLOAD_ARCH}\"")
      elseif(OFFLOAD_TARGET MATCHES "nvptx-none")
        set(OPENMP_OFFLOAD_COMPILE_OPTIONS
            "${OPENMP_OFFLOAD_COMPILE_OPTIONS} -foffload-options=${OFFLOAD_TARGET}=\"-misa=${OFFLOAD_ARCH}\"")
      else()
        message(
          WARNING
            "We don't know how to handle OFFLOAD_ARCH=${OFFLOAD_ARCH} for OFFLOAD_TARGET=${OFFLOAD_TARGET}. Got ignored."
        )
      endif()
    endif()
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -foffload=disable")
  endif()
endif(QMC_OMP)

# Set gnu specific flags (which we always want)
add_definitions(-Drestrict=__restrict__)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -finline-limit=1000 -fstrict-aliasing -funroll-all-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -finline-limit=1000 -fstrict-aliasing -funroll-all-loops")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")

# Suppress compile warnings
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

# treat VLA as error
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror=vla")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wvla")

# set compiler warnings
set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wcomment -Wmisleading-indentation -Wmaybe-uninitialized -Wuninitialized -Wreorder -Wno-unknown-pragmas -Wno-sign-compare"
)

if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.2)
  string(APPEND CMAKE_CXX_FLAGS " -Wsuggest-override")
endif()

# Set extra optimization specific flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

#--------------------------------------
# Special architectural flags
#--------------------------------------
# case arch
#     x86_64: -march
#     powerpc: -mpcu
#     arm: -mpcu
#     default or cray: none
#--------------------------------------
if(CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
  # It's a cray machine. Don't do anything
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  # the case for x86_64
  # check if the user has already specified -march=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # make sure that the user specifies -march= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")

    else() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
      message(
        FATAL_ERROR
          "if -march=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # use -march=native
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  # the case for PowerPC and ARM
  # check if the user has already specified -mcpu=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # make sure that the user specifies -mcpu= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")

    else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
      message(
        FATAL_ERROR
          "if -mcpu=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # use -mcpu=native
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
endif()

# protect buggy libmvec
file(
  WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_glibc.cxx"
  "#include <iostream>\n#if __GLIBC__ == 2 && ( __GLIBC_MINOR__ == 22 || __GLIBC_MINOR__ == 23 )\n#error buggy glibc version\n#endif\n int main() { return 0; }\n"
)
try_compile(PASS_GLIBC ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/src_glibc.cxx
            CMAKE_FLAGS "${CMAKE_CXX_FLAGS}")
if(NOT PASS_GLIBC)
  message(FATAL_ERROR "Your system and GNU compiler are using glibc 2.22 or 2.23 which contains a buggy libmvec."
                      "This results in crashes. Workaround needed. Alternatively upgrade or use another compiler.")
endif()

# Add static flags if necessary
if(QMC_BUILD_STATIC)
  set(CMAKE_CXX_LINK_FLAGS " -static")
endif(QMC_BUILD_STATIC)

# Coverage
if(ENABLE_GCOV)
  set(GCOV_SUPPORTED TRUE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif(ENABLE_GCOV)
