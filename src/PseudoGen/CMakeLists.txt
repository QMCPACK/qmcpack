
PROJECT(ppg)
                                                          
SET(PPSRCS 
  ../Optimize/Minimize.cpp
  ../Numerics/Clebsch_Gordan.cpp
  ../SQD/SphericalPotential/RadialPotential.cpp
  ../SQD/SphericalPotential/SJPseudoPotential.cpp
  PseudoGen.cpp
  PseudoGen.IO.cpp
  PlotPP.cpp
  PPApp.cpp
)

SET(PPSRCS ${PPSRCS})
ADD_EXECUTABLE(ppg ${PPSRCS} SQDFrame.cpp)
TARGET_LINK_LIBRARIES(ppg qmcbase)

