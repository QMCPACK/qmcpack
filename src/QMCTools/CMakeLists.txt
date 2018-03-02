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

SET(QMCTOOLS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  RGFBuilderBase.cpp
  GTO2GridBuilder.cpp
)

SET(MOSRCS
  ../QMCWaveFunctions/OrbitalBuilderBase.cpp
  ../QMCApp/InitMolecularSystem.cpp
  ../QMCApp/ParticleSetPool.cpp
  RGFBuilderBase.cpp
  GTO2GridBuilder.cpp
  QMCGaussianParserBase.cpp
  GaussianFCHKParser.cpp
  GamesXmlParser.cpp
  GamesAsciiParser.cpp
  CasinoParser.cpp
  BParser.cpp
  BMakeFunc.cpp
  GamesFMOParser.cpp
  VSVBParser.cpp
  QPParser.cpp
  PyscfParser.cpp
)

# create libmocommon 
ADD_LIBRARY(mocommon ${MOSRCS})

set(QTOOLS convert4qmc MSDgenerator extract-eshdf-kvectors trace_density)
ADD_EXECUTABLE(getSupercell getSupercell.cpp)

FOREACH(p ${QTOOLS})

  ADD_EXECUTABLE(${p} ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} mocommon qmcbase qmcutil)

  FOREACH(l ${QMC_UTIL_LIBS})
    TARGET_LINK_LIBRARIES(${p} ${l})
  ENDFOREACH(l ${QMC_UTIL_LIBS})

  if(ENABLE_TAU_PROFILE)
    TARGET_LINK_LIBRARIES(${p} tau)
  endif(ENABLE_TAU_PROFILE)

  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)

ENDFOREACH(p ${QTOOLS})

