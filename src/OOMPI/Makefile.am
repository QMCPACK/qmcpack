# -*- makefile -*-
#
# Copyright (c) 2002-2003 Indiana University.  All rights reserved.
# Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
#                         All rights reserved.
# 
# This file is part of the OOMPI software package.  For license
# information, see the LICENSE file in the top level directory of the
# OOMPI source distribution.
#
# $Id$
#

#INCLUDES		= -I$(top_srcdir)/src 

lib_LIBRARIES		= liboompi.a

OO_SRC= \
	Cart_comm.cc \
	Constants.cc \
	Comm.cc \
	Comm_world.cc \
	Datatype.cc \
	Environment.cc \
	Error_table.cc \
	Graph_comm.cc \
	Group.cc \
	Inter_comm.cc \
	Intra_comm.cc \
	Linked_list.cc \
	Message.cc \
	Op.cc \
	Packed.cc \
	Port.cc \
	Request.cc \
	Status.cc \
	User_type.cc \
	Util.cc

OO_INC= Cart_comm.h \
	Comm.h \
	Comm_world.h \
	Constants.h \
	Datatype.h \
	Environment.h \
	Error.h \
	Error_table.h \
	Graph_comm.h \
	Group.h \
	Hidden.h \
	Inter_comm.h \
	Intra_comm.h \
	Linked_list.h \
	Message.h \
	Op.h \
	Packed.h \
	Port.h \
	Request.h \
	Status.h \
	Tag.h \
	User_type.h \
	Util.h \
	oompi.h \
	Wrapper_ptr.cct \
	Wrapper_obj.cct

liboompi_a_SOURCES	= $(OO_SRC)
pkginclude_HEADERS	= $(OO_INC)

