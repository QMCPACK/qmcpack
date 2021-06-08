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

project(qmctools)

set(MOSRCS QMCGaussianParserBase.cpp GaussianFCHKParser.cpp GamesAsciiParser.cpp LCAOHDFParser.cpp)

# create libmocommon
add_library(mocommon ${MOSRCS})
target_link_libraries(mocommon PUBLIC qmcparticle)

set(QTOOLS convert4qmc qmc-extract-eshdf-kvectors)
foreach(p ${QTOOLS})

  add_executable(${p} ${p}.cpp)
  target_link_libraries(${p} mocommon)
  if(USE_OBJECT_TARGET)
    target_link_libraries(${p} qmcparticle qmcutil containers platform_omptarget)
  endif()

  install(TARGETS ${p} RUNTIME DESTINATION bin)

endforeach(p ${QTOOLS})

add_executable(qmc-get-supercell getSupercell.cpp)
install(TARGETS qmc-get-supercell RUNTIME DESTINATION bin)

add_executable(qmc-check-affinity check-affinity.cpp)
if(HAVE_MPI)
  target_link_libraries(qmc-check-affinity MPI::MPI_CXX)
endif()

install(TARGETS qmc-check-affinity RUNTIME DESTINATION bin)

add_executable(convertpw4qmc convertpw4qmc.cpp XmlRep.cpp WriteEshdf.cpp)
target_link_libraries(convertpw4qmc qmcutil Math::FFTW3)

install(TARGETS convertpw4qmc RUNTIME DESTINATION bin)

set(FSSRCS QMCFiniteSize/QMCFiniteSize.cpp QMCFiniteSize/SkParserBase.cpp QMCFiniteSize/SkParserASCII.cpp
           QMCFiniteSize/SkParserScalarDat.cpp QMCFiniteSize/SkParserHDF5.cpp QMCFiniteSize/FSUtilities.cpp)

add_library(fstool ${FSSRCS})
target_link_libraries(fstool PUBLIC qmc qmcparticle qmcwfs)

add_executable(qmcfinitesize qmcfinitesize.cpp)
target_link_libraries(qmcfinitesize fstool)
if(USE_OBJECT_TARGET)
  target_link_libraries(qmcfinitesize qmcparticle qmcutil containers platform_omptarget)
endif()

install(TARGETS qmcfinitesize RUNTIME DESTINATION bin)

add_subdirectory(ppconvert)
