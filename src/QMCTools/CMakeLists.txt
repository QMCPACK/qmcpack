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
  QMCGaussianParserBase.cpp
  GaussianFCHKParser.cpp
  GamesAsciiParser.cpp
  DiracParser.cpp
  CasinoParser.cpp
  BParser.cpp
  BMakeFunc.cpp
  LCAOHDFParser.cpp
)

# create libmocommon 
ADD_LIBRARY(mocommon ${MOSRCS})
TARGET_LINK_LIBRARIES(mocommon PUBLIC qmcparticle)

set(QTOOLS convert4qmc qmc-extract-eshdf-kvectors)
FOREACH(p ${QTOOLS})

  ADD_EXECUTABLE(${p} ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} mocommon)
  IF(USE_OBJECT_TARGET)
    TARGET_LINK_LIBRARIES(${p} qmcparticle qmcutil containers platform_omptarget)
  ENDIF()

  INSTALL(TARGETS ${p}
          RUNTIME DESTINATION bin)

ENDFOREACH(p ${QTOOLS})

ADD_EXECUTABLE(qmc-get-supercell getSupercell.cpp)
INSTALL(TARGETS qmc-get-supercell
        RUNTIME DESTINATION bin)

ADD_EXECUTABLE(qmc-check-affinity check-affinity.cpp)
IF(HAVE_MPI)
  TARGET_LINK_LIBRARIES(qmc-check-affinity MPI::MPI_CXX)
ENDIF()

INSTALL(TARGETS qmc-check-affinity
        RUNTIME DESTINATION bin)

ADD_EXECUTABLE(convertpw4qmc convertpw4qmc.cpp XmlRep.cpp WriteEshdf.cpp)
TARGET_LINK_LIBRARIES(convertpw4qmc qmcutil Math::FFTW3)

INSTALL(TARGETS convertpw4qmc
        RUNTIME DESTINATION bin)

SET(FSSRCS
  QMCFiniteSize/QMCFiniteSize.cpp
  QMCFiniteSize/SkParserBase.cpp
  QMCFiniteSize/SkParserASCII.cpp
  QMCFiniteSize/SkParserScalarDat.cpp
  QMCFiniteSize/SkParserHDF5.cpp
  QMCFiniteSize/FSUtilities.cpp
)

ADD_LIBRARY(fstool ${FSSRCS})
TARGET_LINK_LIBRARIES(fstool PUBLIC qmc qmcparticle qmcwfs)

ADD_EXECUTABLE(qmcfinitesize qmcfinitesize.cpp)
TARGET_LINK_LIBRARIES(qmcfinitesize fstool)
IF(USE_OBJECT_TARGET)
  TARGET_LINK_LIBRARIES(qmcfinitesize qmcparticle qmcutil containers platform_omptarget)
ENDIF()


INSTALL(TARGETS qmcfinitesize
        RUNTIME DESTINATION bin)

#ppconvert
if( BUILD_PPCONVERT )
  SUBDIRS(ppconvert)
endif( BUILD_PPCONVERT )
