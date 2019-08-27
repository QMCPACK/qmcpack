# Check compiler version
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4 )
  MESSAGE(STATUS "Compiler Version ${CMAKE_CXX_COMPILER_VERSION}")
  MESSAGE(FATAL_ERROR "Requires clang 3.4 or higher ")
ENDIF()

# Set the std
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -std=c99")

# Enable OpenMP
IF(QMC_OMP)
  SET(ENABLE_OPENMP 1)
  SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -fopenmp")
  IF(ENABLE_OFFLOAD)
    SET(OFFLOAD_TARGET "nvptx64-nvidia-cuda" CACHE STRING "Offload target architecture")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fopenmp-targets=${OFFLOAD_TARGET}")
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  ENDIF()
ENDIF(QMC_OMP)

# Set clang specific flags (which we always want)
ADD_DEFINITIONS( -Drestrict=__restrict__ )

SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -fomit-frame-pointer -fstrict-aliasing")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -fstrict-aliasing -D__forceinline=inline")

# Suppress compile warnings
SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wno-deprecated -Wno-unused-value")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unused-value")
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.8.0" )
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
ENDIF()

# Set extra optimization specific flags
SET( CMAKE_C_FLAGS_RELEASE     "${CMAKE_C_FLAGS_RELEASE} -ffast-math" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" )
SET( CMAKE_C_FLAGS_RELWITHDEBINFO     "${CMAKE_C_FLAGS_RELWITHDEBINFO} -ffast-math" )
SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ffast-math" )

# Set extra debug flags
SET( CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG} -fno-omit-frame-pointer -fstandalone-debug" )
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fstandalone-debug" )

# Special architectural flags
#--------------------------------------
# case arch
#     x86_64: -march
#     powerpc: -mpcu
#     arm: -mpcu
#     default or cray: none
#--------------------------------------
IF($ENV{CRAYPE_VERSION} MATCHES ".")
  # It's a cray machine. Don't do anything
ELSEIF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  # the case for x86_64
  # check if the user has already specified -march=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # make sure that the user specifies -march= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
    else() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
      MESSAGE(FATAL_ERROR "if -march=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-march=" AND CMAKE_C_FLAGS MATCHES "-march=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
    # use -march=native
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -march=native")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
  endif()  #(CMAKE_CXX_FLAGS MATCHES "-march=" OR CMAKE_C_FLAGS MATCHES "-march=")
ELSEIF(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  # the case for PowerPC and ARM
  # check if the user has already specified -mcpu=XXXX option for cross-compiling.
  if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # make sure that the user specifies -mcpu= for both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS.
    if(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
    else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
      MESSAGE(FATAL_ERROR "if -mcpu=ARCH is specified by the user, it should be added in both CMAKE_CXX_FLAGS and CMAKE_C_FLAGS!")
    endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" AND CMAKE_C_FLAGS MATCHES "-mcpu=")
  else() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
    # use -mcpu=native
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -mcpu=native")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
  endif() #(CMAKE_CXX_FLAGS MATCHES "-mcpu=" OR CMAKE_C_FLAGS MATCHES "-mcpu=")
ENDIF()

# Add static flags if necessary
IF(QMC_BUILD_STATIC)
    SET(CMAKE_CXX_LINK_FLAGS " -static")
ENDIF(QMC_BUILD_STATIC)

# Coverage
IF (ENABLE_GCOV)
  SET(GCOV_COVERAGE TRUE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
ENDIF(ENABLE_GCOV)

SET(XRAY_PROFILE FALSE CACHE BOOL "Use llvm xray profiling")
SET(XRAY_INSTRUCTION_THRESHOLD 200 CACHE STRING "Instruction threshold for xray instrumentation")

IF(XRAY_PROFILE)
  set(XRAY_FLAGS "-fxray-instrument -fxray-instruction-threshold=${XRAY_INSTRUCTION_THRESHOLD}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XRAY_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${XRAY_FLAGS}")
ENDIF(XRAY_PROFILE)

SET(LLVM_SANITIZE_ADDRESS FALSE CACHE BOOL "Use llvm address sanitizer library")
MARK_AS_ADVANCED(LLVM_SANITIZE_ADDRESS)
IF(LLVM_SANITIZE_ADDRESS)
  SET(CMAKE_C_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope ${CMAKE_C_FLAGS_DEBUG}")
  SET(CMAKE_CXX_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope ${CMAKE_CXX_FLAGS_DEBUG}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address -fsanitize-address-use-after-scope")
ENDIF(LLVM_SANITIZE_ADDRESS)

# Don't expect this to be useful unless you have msan instrumented all libraries
SET(LLVM_SANITIZE_MEMORY FALSE CACHE BOOL "Use llvm memory sanitizer library")
MARK_AS_ADVANCED(LLVM_SANITIZE_MEMORY)
IF(LLVM_SANITIZE_MEMORY)
  SET(LLVM_BLACKLIST_SANITIZE_MEMORY "-fsanitize-blacklist=${PROJECT_SOURCE_DIR}/llvm_misc/memory_sanitizer_blacklist.txt")
  SET(CMAKE_C_FLAGS_DEBUG "-fsanitize=memory ${LLVM_BLACKLIST_SANITIZE_MEMORY} ${CMAKE_C_FLAGS_DEBUG}")
  SET(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=memory ${LLVM_BLACKLIST_SANITIZE_MEMORY} ${CMAKE_CXX_FLAGS_DEBUG}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=memory ${LLVM_BLACKLIST_SANITIZE_MEMORY}")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=memory ${LLVM_BLACKLIST_SANITIZE_MEMORY}")
ENDIF(LLVM_SANITIZE_MEMORY)
