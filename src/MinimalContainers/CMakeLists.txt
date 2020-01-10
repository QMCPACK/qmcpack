
SET(MINCONT_SRC ConstantSizeVector.cpp ConstantSizeMatrix.cpp)

ADD_LIBRARY(qmcminimalcontainers ${MINCONT_SRC} )

SUBDIRS(tests)
