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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${qmcpack_BINARY_DIR}/bin)
add_executable(convert4qmc convert4qmc.cpp QMCGaussianParserBase.cpp GaussianFCHKParser.cpp GamesAsciiParser.cpp
  LCAOHDFParser.cpp DiracParser.cpp RMGParser.cpp)

target_link_libraries(convert4qmc PUBLIC qmcparticle)
if(USE_OBJECT_TARGET)
  target_link_libraries(convert4qmc PUBLIC qmcutil)
endif()

add_executable(qmc-extract-eshdf-kvectors qmc-extract-eshdf-kvectors.cpp)
target_link_libraries(qmc-extract-eshdf-kvectors PUBLIC qmcio)

add_executable(qmc-get-supercell getSupercell.cpp)


add_executable(qmc-check-affinity check-affinity.cpp)
if(HAVE_MPI)
  target_link_libraries(qmc-check-affinity MPI::MPI_CXX)
endif()


add_executable(convertpw4qmc convertpw4qmc.cpp XmlRep.cpp WriteEshdf.cpp)
target_link_libraries(convertpw4qmc qmcutil Math::FFTW3)

set(FSSRCS QMCFiniteSize/QMCFiniteSize.cpp QMCFiniteSize/SkParserBase.cpp QMCFiniteSize/SkParserASCII.cpp
           QMCFiniteSize/SkParserScalarDat.cpp QMCFiniteSize/SkParserHDF5.cpp QMCFiniteSize/FSUtilities.cpp)

add_library(fstool ${FSSRCS})
target_link_libraries(fstool PUBLIC qmc qmcparticle qmcwfs)

add_executable(qmcfinitesize qmcfinitesize.cpp)
target_link_libraries(qmcfinitesize fstool)
if(USE_OBJECT_TARGET)
  target_link_libraries(qmcfinitesize qmcparticle qmcutil)
endif()

foreach(EXE_TARGET convert4qmc qmc-extract-eshdf-kvectors qmc-get-supercell qmc-check-affinity convertpw4qmc qmcfinitesize)
  add_test_target_in_output_location(${EXE_TARGET} bin)
  install(TARGETS ${EXE_TARGET} RUNTIME DESTINATION bin)
endforeach()

add_subdirectory(ppconvert)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
