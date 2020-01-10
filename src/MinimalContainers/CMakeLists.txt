
SET(MINCONT_SRC ConstantSizeVector.cpp)

ADD_LIBRARY(qmcminimalcontainers ${MINCONT_SRC} )

SUBDIRS(tests)
