##########################################################################
# Imported from
# http://svn.boost.org/svn/boost/sandbox-branches/boost-cmake/boost_1_34_0
##########################################################################
# Copyright (C) 2007 Douglas Gregor <doug.gregor@gmail.com>
# Copyright (C) 2007 Troy Straszheim
#                                                                        #
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or copy at
#   http://www.boost.org/LICENSE_1_0.txt
##########################################################################
##########################################################################
# Boost Version information
# 1.34.x supports cmake build.
##########################################################################
#set(BOOST_VERSION_MAJOR 1)
#set(BOOST_VERSION_MINOR 33)
#set(BOOST_VERSION_SUBMINOR 0)
#set(BOOST_VERSION "${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}.${BOOST_VERSION_SUBMINOR}")
#SET (BOOST_LIBS $ENV{BOOST_HOME}/libs)

set(Boost_ROOT_DIR ${Boost_INCLUDE_DIR})
set(BOOST_SOURCES "")

if(EXISTS ${Boost_ROOT_DIR}/libs/serialization/src)
  set(BUILD_BOOST_SERIALIZATION TRUE)
endif()

if(EXISTS ${Boost_ROOT_DIR}/libs/filesystems/src)
  set(BUILD_BOOST_FILESYSTEM TRUE)
endif()

if(EXISTS ${Boost_ROOT_DIR}/libs/program_options/src)
  set(BUILD_BOOST_PROGRAM_OPTIONS TRUE)
endif()

set(BUILD_BOOST_MPI FALSE)
if(HAVE_MPI)
  if(EXISTS ${Boost_ROOT_DIR}/libs/mpi/src)
  set(BUILD_BOOST_MPI TRUE)
  endif()
endif(HAVE_MPI)

if(BUILD_BOOST_SERIALIZATION)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/serialization/src")
  set(SOURCES archive_exception.cpp basic_archive.cpp
    basic_iarchive.cpp basic_iserializer.cpp basic_oarchive.cpp
    basic_oserializer.cpp basic_pointer_iserializer.cpp
    basic_pointer_oserializer.cpp
    basic_serializer_map.cpp
    basic_text_iprimitive.cpp
    basic_text_oprimitive.cpp
    basic_text_wiprimitive.cpp
    basic_text_woprimitive.cpp
    basic_xml_archive.cpp
    binary_iarchive.cpp
    binary_oarchive.cpp
    binary_wiarchive.cpp
    binary_woarchive.cpp
    codecvt_null.cpp
    extended_type_info.cpp
    extended_type_info_no_rtti.cpp
    extended_type_info_typeid.cpp
    polymorphic_iarchive.cpp
    polymorphic_oarchive.cpp
    shared_ptr_helper.cpp
    stl_port.cpp
    text_iarchive.cpp
    text_oarchive.cpp
    text_wiarchive.cpp
    text_woarchive.cpp
    utf8_codecvt_facet.cpp
    void_cast.cpp
    xml_archive_exception.cpp
    xml_grammar.cpp
    xml_iarchive.cpp
    xml_oarchive.cpp
    xml_wgrammar.cpp
    xml_wiarchive.cpp
    xml_woarchive.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_SERIALIZATION)

# Boost.Filesystem
if(BUILD_BOOST_FILESYSTEM)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/filesystem/src")
  set(SOURCES
    operations.cpp path.cpp portability.cpp utf8_codecvt_facet.cpp)
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_FILESYSTEM)

if(BUILD_BOOST_MPI)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/mpi/src")
  set(SOURCES
    boost_mpi
    broadcast.cpp
    communicator.cpp
    computation_tree.cpp
    content_oarchive.cpp
    environment.cpp
    exception.cpp
    graph_communicator.cpp
    group.cpp
    intercommunicator.cpp
    mpi_datatype_cache.cpp
    mpi_datatype_oarchive.cpp
    packed_iarchive.cpp
    packed_oarchive.cpp
    packed_skeleton_iarchive.cpp
    packed_skeleton_oarchive.cpp
    point_to_point.cpp
    request.cpp
    text_skeleton_oarchive.cpp
    )
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_MPI)

if(BUILD_BOOST_PROGRAM_OPTIONS)
  set(DIRECTORY "${Boost_ROOT_DIR}/libs/program_options/src")
  set(SOURCES
      cmdline.cpp config_file.cpp options_description.cpp parsers.cpp
      variables_map.cpp value_semantic.cpp positional_options.cpp
      utf8_codecvt_facet.cpp convert.cpp winmain.cpp
      )
  foreach(S ${SOURCES})
    if(EXISTS ${DIRECTORY}/${S})
      set(BOOST_SOURCES ${BOOST_SOURCES} ${DIRECTORY}/${S})
    endif(EXISTS ${DIRECTORY}/${S})
  endforeach(S)
endif(BUILD_BOOST_PROGRAM_OPTIONS)


ADD_LIBRARY(boost ${BOOST_SOURCES})
link_libraries(boost)

#SUBDIRS(serialization)
