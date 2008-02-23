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

SET (BOOST_LIBS $ENV{BOOST_HOME}/libs)

set(SOURCES
    ${BOOST_LIBS}/serialization/src/basic_archive.cpp
    ${BOOST_LIBS}/serialization/src/basic_archive_impl.cpp
    ${BOOST_LIBS}/serialization/src/basic_iarchive.cpp
    ${BOOST_LIBS}/serialization/src/basic_oarchive.cpp
    ${BOOST_LIBS}/serialization/src/basic_iserializer.cpp
    ${BOOST_LIBS}/serialization/src/basic_oserializer.cpp
    ${BOOST_LIBS}/serialization/src/basic_pointer_iserializer.cpp
    ${BOOST_LIBS}/serialization/src/basic_pointer_oserializer.cpp
    ${BOOST_LIBS}/serialization/src/basic_serializer_map.cpp
    ${BOOST_LIBS}/serialization/src/basic_text_iprimitive.cpp
    ${BOOST_LIBS}/serialization/src/basic_text_oprimitive.cpp
    ${BOOST_LIBS}/serialization/src/basic_xml_archive.cpp
    ${BOOST_LIBS}/serialization/src/binary_iarchive.cpp
    ${BOOST_LIBS}/serialization/src/binary_oarchive.cpp
    ${BOOST_LIBS}/serialization/src/extended_type_info.cpp
    ${BOOST_LIBS}/serialization/src/extended_type_info_no_rtti.cpp
    ${BOOST_LIBS}/serialization/src/extended_type_info_typeid.cpp
    ${BOOST_LIBS}/serialization/src/polymorphic_iarchive.cpp
    ${BOOST_LIBS}/serialization/src/polymorphic_oarchive.cpp
    ${BOOST_LIBS}/serialization/src/stl_port.cpp
    ${BOOST_LIBS}/serialization/src/text_iarchive.cpp
    ${BOOST_LIBS}/serialization/src/text_oarchive.cpp
    ${BOOST_LIBS}/serialization/src/void_cast.cpp
    ${BOOST_LIBS}/serialization/src/xml_grammar.cpp
    ${BOOST_LIBS}/serialization/src/xml_iarchive.cpp
    ${BOOST_LIBS}/serialization/src/xml_oarchive.cpp
    )
    
set(WSOURCES
    ${BOOST_LIBS}/serialization/src/codecvt_null.cpp
    ${BOOST_LIBS}/serialization/src/utf8_codecvt_facet.cpp
    ${BOOST_LIBS}/serialization/src/basic_text_wiprimitive.cpp
    ${BOOST_LIBS}/serialization/src/basic_text_woprimitive.cpp
    ${BOOST_LIBS}/serialization/src/binary_wiarchive.cpp
    ${BOOST_LIBS}/serialization/src/binary_woarchive.cpp
    ${BOOST_LIBS}/serialization/src/text_wiarchive.cpp
    ${BOOST_LIBS}/serialization/src/text_woarchive.cpp
    ${BOOST_LIBS}/serialization/src/xml_wgrammar.cpp
    ${BOOST_LIBS}/serialization/src/xml_wiarchive.cpp
    ${BOOST_LIBS}/serialization/src/xml_woarchive.cpp
    )

set(POSRCS
    ${BOOST_LIBS}/program_options/src/winmain.cpp
    ${BOOST_LIBS}/program_options/src/variables_map.cpp
    ${BOOST_LIBS}/program_options/src/value_semantic.cpp
    ${BOOST_LIBS}/program_options/src/utf8_codecvt_facet.cpp
    ${BOOST_LIBS}/program_options/src/positional_options.cpp
    ${BOOST_LIBS}/program_options/src/parsers.cpp
    ${BOOST_LIBS}/program_options/src/options_description.cpp
    ${BOOST_LIBS}/program_options/src/convert.cpp
    ${BOOST_LIBS}/program_options/src/config_file.cpp
    ${BOOST_LIBS}/program_options/src/cmdline.cpp
)

ADD_LIBRARY(boost ${SOURCES} ${WSOURCES} ${POSRCS})
#SUBDIRS(serialization)
