PROJECT(qmcPlusPlus)

#define the paths for library and executable for out-of-soruce compilation
SET (LIBRARY_OUTPUT_PATH ${qmcPlusPlus_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
SET (EXECUTABLE_OUTPUT_PATH ${qmcPlusPlus_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

#list of options to determine libraries. First, try to use libraries available
SET(INTEL_COMPILER 0)
SET(ENABLE_OPENMP 0)
SET(HAVE_MPI 0)
SET(HAVE_OOMPI 0)
SET(HAVE_LIBXML2 0)
SET(HAVE_LIBXMLPP 0)
SET(HAVE_LIBHDF5 0)
SET(HAVE_LIBSPRNG 0)
SET(HAVE_LIBBLITZ 0)
SET(HAVE_LIBGSL 0)

OPTION(QMC_SPRNG "Enable SPRNG random number generation library" OFF)
OPTION(QMC_BLITZ "Enable blitz++ library" ON)

######################################################################
#important settings for QMC simulations to instantiate template classes
#OHMMS_* is to reuse Ohmms class libraries
#OHMMS_DIM =  dimension of the problem
#OHMMS_INDEXTYPE = type of index
#OHMMS_PRECISION  = base precision, float, double etc
#QMC_FASTWALKER  = true, when operation is faster for each walker
#QMC_APPS  =  an application to compile
######################################################################
SET(OHMMS_DIM 3)
SET(OHMMS_INDEXTYPE int)
SET(OHMMS_PRECISION double)
SET(QMC_FASTWALKER 1)
SET(QMCAPPS  MolecuApps)

######################################################################
#experimental settings for optimization
######################################################################
SET(INTEL_COMPILER 0)

SET (PROJECT_CMAKE ${qmcPlusPlus_SOURCE_DIR}/CMake)
######################################################
#COMPILER choose one of the cmake files to customize the compiler options
#If nothing is chosen, default settings by cmake will  be used.
######################################################
#INCLUDE(${PROJECT_CMAKE}/Intel8.cmake)
#INCLUDE(${PROJECT_CMAKE}/Intel7.cmake)
INCLUDE(${PROJECT_CMAKE}/GNUCompilers.cmake)
#INCLUDE(${PROJECT_CMAKE}/IBMCompilers.cmake)
#INCLUDE(${PROJECT_CMAKE}/KCCCompilers.cmake)
#INCLUDE(${PROJECT_CMAKE}/Compilers.cmake)
 
INCLUDE(${CMAKE_ROOT}/Modules/FindThreads.cmake)

#check external libraries: required
INCLUDE(${PROJECT_CMAKE}/FindLibxml2.cmake)
#INCLUDE(${PROJECT_CMAKE}/FindLibxmlpp.cmake)
INCLUDE(${PROJECT_CMAKE}/FindLapack.cmake)
INCLUDE(${PROJECT_CMAKE}/FindHDF5.cmake)
INCLUDE(${PROJECT_CMAKE}/FindGSL.cmake)
INCLUDE(${PROJECT_CMAKE}/FindBlitz++.cmake)
 
IF(FOUND_LIBXML2)
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  LINK_LIBRARIES(${LIBXML2_LIBRARY})
  SET(HAVE_LIBXML2 1)
ENDIF(FOUND_LIBXML2)

#IF(FOUND_LIBXMLPP)
#  INCLUDE_DIRECTORIES(${LIBXMLPP_INCLUDE_DIR})
#  LINK_LIBRARIES(${LIBXMLPP_LIBRARY})
#  SET(HAVE_LIBXMLPP 1)
#ENDIF(FOUND_LIBXMLPP)

IF(FOUND_HDF5)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARY})
ENDIF(FOUND_HDF5)

IF(FOUND_GSL)
  SET(HAVE_LIBGSL 1)
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
  LINK_LIBRARIES(${GSL_LIBRARY})
ENDIF(FOUND_GSL)

#check optional external libraries: enable only when the libraries are found
IF(FOUND_BLITZ)
  INCLUDE_DIRECTORIES(${BLITZ_INCLUDE_DIR})
  SET(HAVE_LIBBLITZ 1)
ENDIF(FOUND_BLITZ)
#IF(QMC_BLITZ)
#ENDIF(QMC_BLITZ)
 
IF(QMC_MPI)
  INCLUDE(${PROJECT_SOURCE_DIR}/CMake/addMPI.cmake)
ENDIF(QMC_MPI)
 
IF(QMC_SPRNG)
  INCLUDE(${PROJECT_SOURCE_DIR}/CMake/addSPRNG.cmake)
ENDIF(QMC_SPRNG)
 
INCLUDE_DIRECTORIES(
  ${qmcPlusPlus_SOURCE_DIR}/src
  ${qmcPlusPlus_BINARY_DIR}/src
)
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

SUBDIRS(src)
