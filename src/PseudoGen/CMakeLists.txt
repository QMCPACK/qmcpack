#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//
#// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////
    
    


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

