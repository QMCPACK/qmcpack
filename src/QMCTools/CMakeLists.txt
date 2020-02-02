#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
#//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
#//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////
    
    



PROJECT(qmctools)

SET(MOSRCS
  ../QMCApp/InitMolecularSystem.cpp
  ../QMCApp/ParticleSetPool.cpp
  QMCGaussianParserBase.cpp
  GaussianFCHKParser.cpp
  GamesXmlParser.cpp
  GamesAsciiParser.cpp
  CasinoParser.cpp
  BParser.cpp
  BMakeFunc.cpp
  GamesFMOParser.cpp
  QPParser.cpp
  LCAOHDFParser.cpp
  convertpw4qmc.cpp 
  WriteEshdf.cpp
  XmlRep.cpp
)

# create libmocommon 
ADD_LIBRARY(mocommon ${MOSRCS})

set(QTOOLS convert4qmc qmc-extract-eshdf-kvectors)
ADD_EXECUTABLE(qmc-get-supercell getSupercell.cpp)
INSTALL(TARGETS qmc-get-supercell
        RUNTIME DESTINATION bin)
ADD_EXECUTABLE(qmc-check-affinity check-affinity.cpp)
INSTALL(TARGETS qmc-check-affinity
        RUNTIME DESTINATION bin)

ADD_EXECUTABLE(convertpw4qmc convertpw4qmc.cpp XmlRep.cpp WriteEshdf.cpp)
FOREACH(l ${QMC_UTIL_LIBS})		       
  TARGET_LINK_LIBRARIES(convertpw4qmc ${l})
ENDFOREACH(l ${QMC_UTIL_LIBS})

INSTALL(TARGETS convertpw4qmc
        RUNTIME DESTINATION bin)


FOREACH(p ${QTOOLS})

  ADD_EXECUTABLE(${p} ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} mocommon qmcbase qmcutil)

  INSTALL(TARGETS ${p}
          RUNTIME DESTINATION bin)

  FOREACH(l ${QMC_UTIL_LIBS})
    TARGET_LINK_LIBRARIES(${p} ${l})
  ENDFOREACH(l ${QMC_UTIL_LIBS})

  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p})
  ENDIF(MPI_LIBRARY)

ENDFOREACH(p ${QTOOLS})

#ppconvert
if( BUILD_PPCONVERT )
  SUBDIRS(ppconvert)
endif( BUILD_PPCONVERT )
