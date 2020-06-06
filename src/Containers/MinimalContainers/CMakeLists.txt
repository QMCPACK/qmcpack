
SET(MINCONT_SRC ConstantSizeMatrix.cpp)

ADD_LIBRARY(qmcminimalcontainers ${MINCONT_SRC} )
TARGET_LINK_LIBRARIES(qmcminimalcontainers PUBLIC containers)
