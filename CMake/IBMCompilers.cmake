if(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le")
  message(STATUS "Power8+ system using xlC/xlc/xlf")

  add_definitions(-Drestrict=__restrict__)

  # Clean up flags

  if(CMAKE_C_FLAGS MATCHES "-qhalt=e")
    set(CMAKE_C_FLAGS "")
  endif()
  if(CMAKE_CXX_FLAGS MATCHES "-qhalt=e")
    set(CMAKE_CXX_FLAGS "")
  endif()

  # Suppress compile warnings
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated -Wno-unused-value")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unused-value")

  # Set extra optimization specific flags
  set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O3")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")

  if(QMC_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qsmp=omp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qsmp=omp")
    if(ENABLE_OFFLOAD)
      set(OPENMP_OFFLOAD_COMPILE_OPTIONS "-qoffload")
    endif(ENABLE_OFFLOAD)
  else(QMC_OMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qnothreaded")
  endif(QMC_OMP)

  # Add static flags if necessary
  if(QMC_BUILD_STATIC)
    set(CMAKE_CXX_LINK_FLAGS " -static")
  endif(QMC_BUILD_STATIC)

endif(CMAKE_SYSTEM_PROCESSOR MATCHES "ppc64le")
