# this module look for lapack/blas and other numerical library support
# it will define the following values
# Since lapack and blas are essential, link_liraries are called.

set(LAPACK_FOUND FALSE)
set(BLAS_FOUND FALSE)
set(MKL_FOUND FALSE)
#
#IF(NOT CMAKE_COMPILER_IS_GNUCXX)

if(${CMAKE_C_COMPILER} MATCHES "icc")
  # Intel composer has everything, 
  if($ENV{MKLROOT} MATCHES "composer")
    include_directories($ENV{MKLROOT}/include)
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
      link_libraries(-L$ENV{MKLROOT}/lib/intel64 -mkl=sequential)
    else()
      link_libraries(-L$ENV{MKLROOT}/lib/ia32 -mkl=sequential)
    endif()
    set(LAPACK_FOUND TRUE)
    set(BLAS_FOUND TRUE)
    set(MKL_FOUND TRUE)
    set(HAVE_MKL TRUE)
    set(HAVE_MKL_VML TRUE)
  else()

    set(mkl_home "")

    if($ENV{MKLROOT} MATCHES "mkl")
      set(mkl_home $ENV{MKLROOT})
    else()
      if($ENV{MKL_HOME} MATCHES "mkl")
        set(mkl_home $ENV{MKL_HOME})
      endif()
    endif()

    if(mkl_home MATCHES "mkl")

      #default MKL libraries 
      STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" MKL_VERSION ${mkl_home})

      if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
        if(${MKL_VERSION} MATCHES "10\\.3\\.[0-4]")
          link_libraries(-L${mkl_home}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm)
        else()
          if(${MKL_VERSION} MATCHES "10\\.[0-2]\\.[0-4]")
            link_libraries(-L${mkl_home}/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm)
          else()
            link_libraries(-L${mkl_home}/lib/em64t -lmkl_lapack -lmkl -lguide)
          endif()
        endif()
      endif()

      if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386")
        if(${MKL_VERSION} MATCHES "10\\.[0-3]\\.[0-4]")
          link_libraries(-L${mkl_home}/lib/ia32 -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm)
        else()
          link_libraries(-L${mkl_home}/lib/ia32 -lmkl_lapack -lmkl -lguide)
        endif()
      endif()

      if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")
        link_libraries(-L${mkl_home}/lib/64 -lmkl_lapack -lmkl -lguide)
      endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ia64")

      set(LAPACK_FOUND TRUE)
      set(BLAS_FOUND TRUE)
      set(MKL_FOUND TRUE)

      FIND_PATH(MKL_INCLUDE_DIR mkl.h ${mkl_home}/include)
      if(MKL_INCLUDE_DIR)
        MESSAGE(STATUS "Header files of MKL libraries are found at " ${MKL_INCLUDE_DIR})
        INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
        set(HAVE_MKL TRUE)
        find_file(mkl_vml_file mkl_vml.h ${mkl_home}/include)
        if(mkl_vml_file)
          set(HAVE_MKL_VML TRUE)
        endif(mkl_vml_file)
      endif(MKL_INCLUDE_DIR)

    endif()
  endif()
else()

  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    SET(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -framework Accelerate")
    SET(LAPACK_LIBRARY_INIT 1 CACHE BOOL "use Mac Framework")
    SET(MAC_VECLIB 1 CACHE BOOL "use Mac Framework")
    MESSAGE(STATUS "Using Framework on Darwin.")
    set(LAPACK_FOUND TRUE)
    set(BLAS_FOUND TRUE)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  if(${CMAKE_SYSTEM_NAME} MATCHES "AIX")
    link_libraries(-lessl -lmass -lmassv)
    set(BLAS_FOUND TRUE)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "AIX")

  if($ENV{LAPACK} MATCHES "lapack")
    link_libraries($ENV{LAPACK})
    set(LAPACK_FOUND TRUE)
  endif($ENV{LAPACK} MATCHES "lapack")

  IF($ENV{ATLAS} MATCHES "atlas")
    # COULD SEARCH THESE but..... 
    set(atlas_libs "lapack;f77blas;cblas;atlas")
    set(LAPACK_FOUND TRUE)
    set(BLAS_FOUND TRUE)
    link_libraries($ENV{ATLAS})
  endif($ENV{ATLAS} MATCHES "atlas")

endif()

if(LAPACK_FOUND AND BLAS_FOUND)
  MESSAGE(STATUS "LAPACK and BLAS libraries are linked to all the applications")
else(LAPACK_FOUND AND BLAS_FOUND)
  MESSAGE(STATUS "Failed to link LAPACK, BLAS, ATLAS libraries with environments. Going to search standard paths.")
  find_library(LAPACK_LIBRARIES lapack)
  find_library(BLAS_LIBRARIES blas)
  if(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
    MESSAGE(STATUS "LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
    MESSAGE(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
    link_libraries(${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    set(LAPACK_FOUND TRUE)
    set(BLAS_FOUND TRUE)
  endif(LAPACK_LIBRARIES AND BLAS_LIBRARIES)
endif(LAPACK_FOUND AND BLAS_FOUND)

#MARK_AS_ADVANCED(
#  LAPACK_LIBRARIES 
#  BLAS_LIBRARIES 
#  )
#IF(USE_SCALAPACK)
#  SET(PNPATHS 
#    ${MKL_PATHS}
#    ${BLACS_HOME}/lib
#    ${SCALAPACK_HOME}/lib
#    /usr/lib
#    /opt/lib
#    /usr/local/lib
#    /sw/lib
#    )
#
#  IF(INTEL_MKL)
#    FIND_LIBRARY(BLACSLIB mkl_blacs_${PLAT}_lp${QMC_BITS} PATHS  ${PNPATHS})
#    FIND_LIBRARY(SCALAPACKLIB mkl_scalapack PATHS  ${PNPATHS})
#  ENDIF(INTEL_MKL)
#
#  IF(NOT SCALAPACKLIB)
#    FIND_LIBRARY(BLACSLIB blacs_MPI-${PLAT}-{BLACSDBGLVL} PATHS  ${PNPATHS})
#    FIND_LIBRARY(BLACSCINIT blacsCinit_MPI-${PLAT}-{BLACSDBGLVL} PATHS  ${PNPATHS})
#    FIND_LIBRARY(SCALAPACKLIB scalapack PATHS  ${PNPATHS})
#  ENDIF(NOT SCALAPACKLIB)
#
#  IF(BLACSLIB AND SCALAPACKLIB)
#    SET(FOUND_SCALAPACK 1 CACHE BOOL "Found scalapack library")
#  ELSE(BLACSLIB AND SCALAPACKLIB)
#    SET(FOUND_SCALAPACK 0 CACHE BOOL "Mising scalapack library")
#  ENDIF(BLACSLIB AND SCALAPACKLIB)
#
#  MARK_AS_ADVANCED(
#    BLACSCINIT
#    BLACSLIB
#    SCALAPACKLIB
#    FOUND_SCALAPACK
#    )
#ENDIF(USE_SCALAPACK)

