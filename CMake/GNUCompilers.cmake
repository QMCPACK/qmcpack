#GNU compilers
IF(CMAKE_COMPILER_IS_GNUCXX) 

  exec_program(${CMAKE_C_COMPILER} ARGS "-dumpversion" OUTPUT_VARIABLE _gcc_version_info)
  string(REGEX REPLACE "^([0-9]+).*$"                   "\\1" GCC_MAJOR ${_gcc_version_info})
  string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*$"          "\\1" GCC_MINOR ${_gcc_version_info})
  string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" GCC_PATCH ${_gcc_version_info})


  if(${GCC_MAJOR} LESS 4)
    MESSAGE(FATAL_ERROR "Require gcc 4.4 or higher ")
  endif()
  if(${GCC_MINOR} LESS 4)
    MESSAGE(FATAL_ERROR "Require gcc 4.4 or higher ")
  endif()

  SET(ENABLE_OPENMP 1)

  ADD_DEFINITIONS(-Drestrict=__restrict__ -DADD_ -DINLINE_ALL=inline)
#  SET(CMAKE_CXX_FLAGS "-O3 -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math -Wno-deprecated ")
#  SET(CMAKE_CXX_FLAGS "-g -O3 -ftemplate-depth-60 -Drestrict=__restrict__ -funroll-all-loops   -finline-limit=1000 -Wno-deprecated ")
  set(GNU_FLAGS "-malign-double -fomit-frame-pointer -ffast-math -fopenmp -O3 -finline-limit=1000 -fstrict-aliasing -funroll-all-loops -Wno-deprecated ")
  SET(CMAKE_CXX_FLAGS "${GNU_FLAGS}") # -ftemplate-depth=60")
  SET(CMAKE_C_FLAGS "${GNU_FLAGS} -std=c99")

  IF(HAVE_POSIX_MEMALIGN)
  SET(CMAKE_TRY_GNU_CC_FLAGS "-mmmx")
  CHECK_C_COMPILER_FLAG(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
  IF(GNU_CC_FLAGS)
    SET(HAVE_MMX 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmmx")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmmx")
  ENDIF(GNU_CC_FLAGS)

  SET(CMAKE_TRY_GNU_CC_FLAGS "-msse")
  CHECK_C_COMPILER_FLAG(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
  IF(GNU_CC_FLAGS)
    SET(HAVE_SSE 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
  ENDIF(GNU_CC_FLAGS)

  SET(CMAKE_TRY_GNU_CXX_FLAGS "-msse2")
  CHECK_C_COMPILER_FLAG(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
  IF(GNU_CC_FLAGS)
    SET(HAVE_SSE2 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
  ENDIF(GNU_CC_FLAGS)

  SET(CMAKE_TRY_GNU_CC_FLAGS "-msse3")
  CHECK_C_COMPILER_FLAG(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
  IF(GNU_CC_FLAGS)
    SET(HAVE_SSE3 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3")
  ENDIF(GNU_CC_FLAGS)

  SET(CMAKE_TRY_GNU_CC_FLAGS "-msse4.1")
  CHECK_C_COMPILER_FLAG(${CMAKE_TRY_GNU_CC_FLAGS} GNU_CC_FLAGS)
  IF(GNU_CC_FLAGS)
    SET(HAVE_SSE41 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
  ENDIF(GNU_CC_FLAGS)
  ENDIF(HAVE_POSIX_MEMALIGN)



  #  SET(CMAKE_CXX_FLAGS "-O6 -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -funroll-all-loops   -finline-limit=1000 -ffast-math -Wno-deprecated -pg")
  #  SET(CMAKE_CXX_FLAGS "-g -ftemplate-depth-60 -Drestrict=__restrict__ -fstrict-aliasing -Wno-deprecated")

  IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -faltivec -framework Accelerate -bind_at_load")
    SET(F77 xlf)
    SET(F77FLAGS -O3)
  ELSE(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  #  SET(FORTRAN_LIBS "-lg2c")
    SET(F77 g77)
    SET(F77FLAGS  -funroll-loops -O3)
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  IF(QMC_BUILD_STATIC)
    SET(CMAKE_CXX_LINK_FLAGS " -static")
  ENDIF(QMC_BUILD_STATIC)

  SET(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  SET(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${CMAKE_C_FLAGS}")

ENDIF(CMAKE_COMPILER_IS_GNUCXX) 

#IF(APPLE)
#  INCLUDE_DIRECTORIES(/sw/include)
#ENDIF(APPLE)
