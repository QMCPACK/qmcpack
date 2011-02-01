PROJECT(benchmark)

#ADD_EXECUTABLE( fft1d  fft1d.cpp)
#TARGET_LINK_LIBRARIES(fft1d qmcutil)
#
#ADD_EXECUTABLE( fft1d_bench  fft1d_bench.cpp)
#TARGET_LINK_LIBRARIES(fft1d_bench qmcutil)

ADD_EXECUTABLE( numerics  numerics.cpp)
TARGET_LINK_LIBRARIES(numerics qmcutil)
