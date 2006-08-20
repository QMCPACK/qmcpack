PROJECT(sandbox)

LINK_LIBRARIES(${FORTRAN_LIBS} 
${LAPACK_LIBRARY} 
${BLAS_LIBRARY} 
)

IF(INTEL_MKL)
  LINK_LIBRARIES(${INTEL_GUIDE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
ELSE(INTEL_MKL)
  LINK_LIBRARIES(${CMAKE_THREAD_LIBS_INIT})
ENDIF(INTEL_MKL)

ADD_EXECUTABLE(mpi_group group.cpp)
TARGET_LINK_LIBRARIES(mpi_group qmc qmcbase qmcutil)

#FIND_LIBRARY(BOOST_SERIALIZATION_LIBRARY  boost_serialization-il-1_33
#/u/ncsa/jnkim/share/boost/bin/boost/libs/serialization/build/libboost_serialization.so/intel-linux/release/shared-linkable-true/
#)

ADD_EXECUTABLE(serialization serialization.cpp)
TARGET_LINK_LIBRARIES(serialization boost)

#ADD_EXECUTABLE(reduce_test reduce_performance_test.cpp)
#TARGET_LINK_LIBRARIES(reduce_test boost)
