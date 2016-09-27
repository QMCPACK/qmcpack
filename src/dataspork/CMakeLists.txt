#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
#//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
#//
#// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
#//////////////////////////////////////////////////////////////////////////////////////
    
    

SET(SRCS 
ScalarObservable.cpp 
ScalarDataSet.cpp 
ScalarDataSetManager.cpp 
DataSporkAnalysis.cpp
)

ADD_EXECUTABLE(datasporkpp ${SRCS} datasporkpp.cpp)
TARGET_LINK_LIBRARIES(datasporkpp boost qmcutil)
