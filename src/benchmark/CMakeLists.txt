PROJECT(benchmark)

ADD_EXECUTABLE( fft1d  fft1d.cpp)
TARGET_LINK_LIBRARIES(fft1d qmcbase qmcutil)
