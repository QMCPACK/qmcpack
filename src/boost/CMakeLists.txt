##########################################################################
# boost.mpi
# Copyright © 2005-2007 Douglas Gregor, Matthias Troyer, Trustees of Indiana University
# Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at htt
##########################################################################
##########################################################################
# Boost Version information
# imported 1.60.0 
#  boost/libs/serialization/src -> boost/libs/serialization
#  boost/libs/mpi/src -> boost/libs/mpi
##########################################################################

set(SOURCES 
  basic_serializer_map.cpp
  binary_oarchive.cpp
  extended_type_info.cpp
  extended_type_info_typeid.cpp
  polymorphic_iarchive.cpp
  polymorphic_oarchive.cpp
  stl_port.cpp
  void_cast.cpp
  basic_iarchive.cpp 
  basic_oarchive.cpp
  basic_oserializer.cpp 
  basic_iserializer.cpp 
  basic_archive.cpp
  archive_exception.cpp 
  )

foreach(S ${SOURCES})
  set(BOOST_SOURCES ${BOOST_SOURCES} libs/serialization/${S})
endforeach(S)

set(SOURCES
  timer.cpp
  text_skeleton_oarchive.cpp
  request.cpp
  point_to_point.cpp
  packed_skeleton_oarchive.cpp
  packed_skeleton_iarchive.cpp
  packed_oarchive.cpp
  packed_iarchive.cpp
  mpi_datatype_oarchive.cpp
  mpi_datatype_cache.cpp
  intercommunicator.cpp
  group.cpp
  graph_communicator.cpp
  environment.cpp
  content_oarchive.cpp
  computation_tree.cpp
  communicator.cpp
  broadcast.cpp
  exception.cpp
)

foreach(S ${SOURCES})
  set(BOOST_SOURCES ${BOOST_SOURCES} libs/mpi/${S})
endforeach(S)

ADD_LIBRARY(boost ${BOOST_SOURCES})
