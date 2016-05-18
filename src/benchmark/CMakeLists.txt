PROJECT(benchmark)

SET(BENCH numerics)

#ADD_EXECUTABLE( fft1d  fft1d.cpp)
#TARGET_LINK_LIBRARIES(fft1d qmcutil)
#
#ADD_EXECUTABLE( fft1d_bench  fft1d_bench.cpp)
#TARGET_LINK_LIBRARIES(fft1d_bench qmcutil)

FOREACH(p ${BENCH})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmcutil)
  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)
ENDFOREACH(p ${BENCH})
