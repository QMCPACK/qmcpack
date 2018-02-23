PROJECT(ppconvert)


#####################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

######################################################################
# Standard test
######################################################################
SET(QMC_CMAKE "/Users/epd/codes/qmcpack/CMake")
SET(CMAKE_MODULE_PATH ${QMC_CMAKE} ${CMAKE_MODULE_PATH})

######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
#If the automatic detection does not work, comment out everything 
#upto COMPILER MANUAL SELECTION and use one of the customized cmake file.
######################################################


INCLUDE(${CMAKE_ROOT}/Modules/CheckCCompilerFlag.cmake)
INCLUDE(${CMAKE_ROOT}/Modules/CheckCXXCompilerFlag.cmake)

  ######################################################
  # Try to identify CPU identity
  ######################################################
  SET(CPU_IDENTITY "generic")
  #INCLUDE(${QMC_CMAKE}/CheckProcessorID.cmake)

  #------------------------------------
  # Identify the compiler -- I don't get why this COMPILER variable is added --PD
  #------------------------------------
  IF ( CMAKE_C_COMPILER_WORKS OR CMAKE_C_COMPILER_WORKS )
    IF( (CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR ( CMAKE_C_COMPILER_ID MATCHES "GNU") )
      SET( COMPILER GNU )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "XL") OR ( CMAKE_C_COMPILER_ID MATCHES "XL") )
      SET( COMPILER IBM )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "Intel") OR (CMAKE_C_COMPILER_ID MATCHES "Intel") )
      SET( COMPILER Intel )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "PGI") OR (CMAKE_C_COMPILER_ID MATCHES "PGI") )
      SET( COMPILER PGI )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "Cray") OR (CMAKE_C_COMPILER_ID MATCHES "Cray") )
      SET( COMPILER Cray )
    ELSEIF( (CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR (CMAKE_C_COMPILER_ID MATCHES "Clang") )
      SET( COMPILER Clang )
    ELSE()
        SET(USING_DEFAULT TRUE)
        MESSAGE("${CMAKE_C_COMPILER_ID}")
        MESSAGE(WARNING "Unknown C/C++ compiler, default flags will be used")
    ENDIF()
  ELSE()
    MESSAGE(WARNING "No compiler identified")
  ENDIF()
  MESSAGE(STATUS "Compiler: ${COMPILER}")

  #------------------------------------
  # Include compiler-specific cmake file
  #------------------------------------
  IF( ${COMPILER} MATCHES "IBM" )
    INCLUDE(${QMC_CMAKE}/IBMCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "Intel" )
    INCLUDE(${QMC_CMAKE}/IntelCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "GNU" )
    INCLUDE(${QMC_CMAKE}/GNUCompilers.cmake)
  ELSEIF( ${COMPILER} MATCHES "Clang" )
    INCLUDE(${QMC_CMAKE}/ClangCompilers.cmake)
  ELSE()
    MESSAGE(WARNING "No default file for compiler (${COMPILER})")
  ENDIF()

# Require C++11 while leaving possibility to fall back by modifying this file
SET(CXX11_NEEDED TRUE)
IF(CXX11_NEEDED)
  SET(CXX11_FLAG "-std=c++11")

  FILE( WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/std11.cxx"
    "int main() { return 0; }\n" )
  try_compile(CXX_COMPILER_SUPPORT_CXX11 ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/std11.cxx
  CMAKE_FLAGS
  OUTPUT_VARIABLE STD11_OUT)

  #check if the CXX compiler supports -std=c++11 option
  # Force the flag on Cray with Intel compiler, because the Cray wrapper
  #  prints an warning that interferes with the flag detection code
  #  with older versions of CMake.
  IF($ENV{CRAYPE_VERSION} MATCHES ".")
    IF( ${COMPILER} MATCHES "Intel" AND NOT CXX_COMPILER_SUPPORT_CXX11)
      SET(CXX_COMPILER_SUPPORT_CXX11 TRUE)
      MESSAGE(STATUS "Forcing C++11 support on Cray with Intel")
    ENDIF()
  ENDIF()

  IF (CXX_COMPILER_SUPPORT_CXX11)
    # avoid repeated -std=c++11 flag
    STRING(REPLACE "++" "\\+\\+" CXX11_FLAG_MATCH ${CXX11_FLAG})
    IF(NOT CMAKE_CXX_FLAGS MATCHES ${CXX11_FLAG_MATCH})
      SET (CMAKE_CXX_FLAGS "${CXX11_FLAG} ${CMAKE_CXX_FLAGS}")
    ENDIF()
  ELSE()
    MESSAGE("CXX_COMPILER_SUPPORT_CXX11 test failed: ${STD11_OUT}")
    MESSAGE(FATAL_ERROR "A compiler supporting C++11 is required. Use a newer C++ compiler.")
  ENDIF()
ENDIF(CXX11_NEEDED)

  #-------------------------------------------------------------------
  # check lapack/blas
  #-------------------------------------------------------------------
  # Quick check for MKL be intel compiler or have set MKL_ROOT
  INCLUDE(${QMC_CMAKE}/FindMKL.cmake)
  IF ( HAVE_MKL )
    # We found MKL, add the link library
    SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${MKL_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAGS}")
    SET(LAPACK_FOUND 1)
    SET(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    SET(LAPACK_LINKER_FLAGS "${MKL_LINKER_FLAGS}")
    #SET(LAPACK_LINKER_FLAGS_$<CONFIG> "${MKL_LINKER_FLAGS}")
  ELSE()	 
    # Search for Blas/Lapack
    find_package(LAPACK)
    # When cmake finds and sets up the MKL libraries with a non-Intel compiler, set the flag
    # This should be tidied up when FindMKL.cmake supports non-Intel compilers.
    IF("${LAPACK_LIBRARIES}" MATCHES "mkl")
      MESSAGE(STATUS "MKL found via LAPACK/BLAS")
      SET( MKL_FOUND 1 )
      SET( HAVE_MKL 1 )
    ENDIF()
  ENDIF()
  IF(LAPACK_FOUND)
    MESSAGE(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
    MESSAGE(STATUS "LAPACK linker flags: ${LAPACK_LINKER_FLAGS}")
    SET(LAPACK_LIBRARY ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
  ELSE()
    MESSAGE(FATAL_ERROR "Could not find required libraries LAPACK &/or BLAS") 
  ENDIF()

MESSAGE(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")  
MESSAGE(STATUS "CXX compiler flags: ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "CC  compiler flags: ${CMAKE_C_FLAGS}")
MESSAGE(STATUS "CXX exe linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

INCLUDE_DIRECTORIES( "/usr/local/opt/llvm/include" ${CMAKE_SOURCE_DIR}/src )
###################################################

###################################################
# Done searching libraries. Start building ppconvert
###################################################

MESSAGE(STATUS "Ready to build ppconvert")

SUBDIRS(src)
