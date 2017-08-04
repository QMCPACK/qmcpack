# Check compiler version
IF ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.0 )
  MESSAGE(STATUS "Compiler Version ${CMAKE_CXX_COMPILER_VERSION}")
  MESSAGE(FATAL_ERROR "Require clang 3.0 or higher ")
ENDIF()

# Enable OpenMP
SET(ENABLE_OPENMP 1)
IF ( ENABLE_OPENMP )
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fopenmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF()

# Set the std
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=c99")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Set clang specfic flags (which we always want)
ADD_DEFINITIONS( -Drestrict=__restrict__ )
ADD_DEFINITIONS( -DADD_ )
ADD_DEFINITIONS( -DINLINE_ALL=inline )
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -fomit-frame-pointer -fstrict-aliasing")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fomit-frame-pointer -fstrict-aliasing -D__forceinline=inline")
SET( HAVE_POSIX_MEMALIGN 0 )    # Clang doesn't support -malign-double

# Suppress compile warnings
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wno-deprecated -Wno-unused-value")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-unused-value -Wno-undefined-var-template")

# Set extra optimization specific flags
SET( CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -ffast-math" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math" )
SET( CMAKE_C_FLAGS_RELWITHDEBINFO   "${CMAKE_C_FLAGS_RELWITHDEBINFO}   -ffast-math" )
SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -ffast-math" )

# Enable mmx/sse instructions if posix_memalign exists
IF(HAVE_POSIX_MEMALIGN)

    # Check for mmx flags
    SET(CMAKE_TRY_CC_FLAGS "-mmmx")
    CHECK_C_COMPILER_FLAG(${CMAKE_TRY_CC_FLAGS} CC_FLAGS)
    IF(CC_FLAGS)
        SET(HAVE_MMX 1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmmx")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmmx")
    ENDIF(CC_FLAGS)

    # Check for msse flags
    SET(CMAKE_TRY_CC_FLAGS "-msse")
    CHECK_C_COMPILER_FLAG(${CMAKE_TRY_CC_FLAGS} CC_FLAGS)
    IF(CC_FLAGS)
        SET(HAVE_SSE 1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
    ENDIF(CC_FLAGS)

    # Check for msse2 flags
    SET(CMAKE_TRY_GNU_CXX_FLAGS "-msse2")
    CHECK_C_COMPILER_FLAG(${CMAKE_TRY_CC_FLAGS} CC_FLAGS)
    IF(CC_FLAGS)
        SET(HAVE_SSE2 1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    ENDIF(CC_FLAGS)

    # Check for msse3 flags
    SET(CMAKE_TRY_CC_FLAGS "-msse3")
    CHECK_C_COMPILER_FLAG(${CMAKE_TRY_CC_FLAGS} CC_FLAGS)
    IF(CC_FLAGS)
        SET(HAVE_SSE3 1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse3")
    ENDIF(CC_FLAGS)

    # Check for msse4.1 flags
    SET(CMAKE_TRY_CC_FLAGS "-msse4.1")
    CHECK_C_COMPILER_FLAG(${CMAKE_TRY_CC_FLAGS} CC_FLAGS)
    IF(CC_FLAGS)
        SET(HAVE_SSE41 1)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
    ENDIF(CC_FLAGS)

ENDIF(HAVE_POSIX_MEMALIGN)

# Add static flags if necessary
IF(QMC_BUILD_STATIC)
    SET(CMAKE_CXX_LINK_FLAGS " -static")
ENDIF(QMC_BUILD_STATIC)

# Add enviornmental flags
SET(CMAKE_CXX_FLAGS "$ENV{CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
SET(CMAKE_C_FLAGS "$ENV{CC_FLAGS} ${CMAKE_C_FLAGS}")

# Coverage
IF (ENABLE_GCOV)
  SET(GCOV_COVERAGE TRUE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
ENDIF(ENABLE_GCOV)

