#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2020 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(IO_XML_SRC Libxml2Doc.cpp XMLParsingString.cpp)

add_library(qmcio_xml ${IO_XML_SRC})

#TARGET_INCLUDE_DIRECTORIES(qmcio_xml PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(qmcio_xml PUBLIC containers LibXml2::LibXml2)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
