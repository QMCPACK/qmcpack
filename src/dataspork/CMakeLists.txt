
SET(SRCS 
ScalarObservable.cpp 
ScalarDataSet.cpp 
ScalarDataSetManager.cpp 
DataSporkAnalysis.cpp
)

ADD_EXECUTABLE(datasporkpp ${SRCS} datasporkpp.cpp)
TARGET_LINK_LIBRARIES(datasporkpp boost qmcutil)
