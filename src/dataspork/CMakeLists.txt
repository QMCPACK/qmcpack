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
    
    

SET(SRCS 
ScalarObservable.cpp 
ScalarDataSet.cpp 
ScalarDataSetManager.cpp 
DataSporkAnalysis.cpp
)

ADD_EXECUTABLE(datasporkpp ${SRCS} datasporkpp.cpp)
TARGET_LINK_LIBRARIES(datasporkpp boost qmcutil)
