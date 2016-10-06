#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////


PROJECT(sandbox)

SET(SANDBOX)

#IF(INTEL_MKL)
#  LINK_LIBRARIES(${INTEL_GUIDE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
#ELSE(INTEL_MKL)
#  LINK_LIBRARIES(${CMAKE_THREAD_LIBS_INIT})
#ENDIF(INTEL_MKL)

#IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#  LINK_LIBRARIES(${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${FORTRAN_LIBS})
#ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#
#IF(HAVE_MPI)
#  ADD_EXECUTABLE(mpi_group group.cpp)
#  TARGET_LINK_LIBRARIES(mpi_group qmc qmcbase qmcutil)
#
#  ADD_EXECUTABLE( observable_helper_test  observable_helper_test.cpp)
#  TARGET_LINK_LIBRARIES( observable_helper_test qmcbase qmcutil)
#ENDIF(HAVE_MPI)

set(MYTEST lattice)
FOREACH(p ${MYTEST})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmcbase qmcutil)

  FOREACH(l ${QMC_UTIL_LIBS})
    TARGET_LINK_LIBRARIES(${p} ${l})
  ENDFOREACH(l ${QMC_UTIL_LIBS})

  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)
ENDFOREACH(p ${MYTEST})


IF(HAVE_EINSPLINE)

  SET(ESTEST einspline_bench einspline_smp einspline_validation)

  FOREACH(p ${ESTEST})
    ADD_EXECUTABLE( ${p}  ${p}.cpp)
    TARGET_LINK_LIBRARIES(${p} qmcbase qmcutil)
    IF(HAVE_EINSPLINE_EXT)
      TARGET_LINK_LIBRARIES(${p} ${EINSPLINE_LIBRARIES})
    else()
      TARGET_LINK_LIBRARIES(${p} einspline)
    endif()
    IF(MPI_LIBRARY)
      TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
    ENDIF(MPI_LIBRARY)

  ENDFOREACH(p ${ESTEST})

ENDIF(HAVE_EINSPLINE)

#ADD_EXECUTABLE(gnubug gnubug.cpp)
#TARGET_LINK_LIBRARIES(gnubug qmcbase qmcutil)
#FIND_LIBRARY(BOOST_SERIALIZATION_LIBRARY  boost_serialization-il-1_33
#/u/ncsa/jnkim/share/boost/bin/boost/libs/serialization/build/libboost_serialization.so/intel-linux/release/shared-linkable-true/
#)

#ADD_EXECUTABLE(serialization serialization.cpp)
#TARGET_LINK_LIBRARIES(serialization boost)

#ADD_EXECUTABLE(cubicspline cubicspline.cpp)
#TARGET_LINK_LIBRARIES(cubicspline qmcbase qmcutil)
#
#ADD_EXECUTABLE(cubicspline1 cubicspline1.cpp)
#TARGET_LINK_LIBRARIES(cubicspline1 qmcbase qmcutil)
#
#ADD_EXECUTABLE(bspline bspline.cpp)
#TARGET_LINK_LIBRARIES(bspline qmc qmcbase)
#
#ADD_EXECUTABLE(tribsplineComplex tribsplineComplex.cpp)
#TARGET_LINK_LIBRARIES(tribsplineComplex qmc qmcbase)
#ADD_EXECUTABLE(tribsplineReal tribsplineReal.cpp)
#TARGET_LINK_LIBRARIES(tribsplineReal qmc qmcbase)


#INCLUDE(${PROJECT_CMAKE}/FindBlitz++.cmake)
#INCLUDE_DIRECTORIES(${BLITZ_INCLUDE_DIR})
#ADD_EXECUTABLE(bspline_blitz bspline_blitz.cpp)
#TARGET_LINK_LIBRARIES(bspline_blitz qmc qmcbase)

#ADD_EXECUTABLE(reduce_test reduce_performance_test.cpp)
#TARGET_LINK_LIBRARIES(reduce_test boost)
