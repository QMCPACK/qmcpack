#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2021 QMCPACK developers.
#//
#// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#//
#// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//////////////////////////////////////////////////////////////////////////////////////

set(IO_XML_SRC Libxml2Doc.cpp XMLParsingString.cpp ParameterSet.cpp libxmldefs.cpp)

add_library(qmcio_xml ${IO_XML_SRC})

# headers used both in and out of the module include based on src/io being the root for
# io header inclusion
target_include_directories(qmcio_xml PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(qmcio_xml PUBLIC cxx_helpers containers LibXml2::LibXml2)

if(BUILD_UNIT_TESTS)
  add_subdirectory(tests)
endif()
