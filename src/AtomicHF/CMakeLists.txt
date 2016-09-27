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
    
    

PROJECT(hf)

#CONFIGURE_FILE(${hf_SOURCE_DIR}/../Configuration.h.in
#              ${hf_SOURCE_DIR}/../Configuration.h)


SET(HFSRCS 
  ../Utilities/OhmmsInform.cpp 
  ../OhmmsApp/ProjectData.cpp 
  ../Platforms/sysutil.cpp
  ../AtomicHF/Clebsch_Gordan.cpp
  ../AtomicHF/RadialPotential.cpp
  ../AtomicHF/HartreeFock.cpp
)

LINK_LIBRARIES(${FORTRAN_LIBS})
ADD_EXECUTABLE(hf ${HFSRCS})
TARGET_LINK_LIBRARIES(hf qmcbase)
