// -*- c++ -*-
//
// Copyright (c) 2002-2003 Indiana University.  All rights reserved.
// Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
//                         All rights reserved.
//
// This file is part of the OOMPI software package.  For license
// information, see the LICENSE file in the top level directory of the
// OOMPI source distribution.
//
// $Id$
//
// OOMPI Class library
//

#ifndef _OOMPI_H_
#define _OOMPI_H_

//
// MPI include file
//

#include <mpi.h>


//
// Forward references
//

class OOMPI_Cart_comm;
class OOMPI_Comm;
class OOMPI_Comm_world;
class OOMPI_Datatype;
class OOMPI_Environment;
class OOMPI_Graph_comm;
class OOMPI_Group;
class OOMPI_Inter_comm;
class OOMPI_Intra_comm;
class OOMPI_Op;
class OOMPI_Packed;
class OOMPI_Message;
class OOMPI_Array_message;
class OOMPI_Port;
class OOMPI_Request;
class OOMPI_Status;
class OOMPI_Tag;
class OOMPI_User_type;


//
// All the individual include files
//


#include <oompi-config.h>
#include "Error.h"
#include "Constants.h"
#include "Tag.h"

#include "Comm.h"
#include "Intra_comm.h"
#include "Inter_comm.h"
#include "Cart_comm.h"
#include "Graph_comm.h"
#include "Comm_world.h"
#include "Port.h"
#include "Group.h"
#include "Request.h"
#include "Status.h"
#include "Packed.h"
#include "Message.h"
#include "User_type.h"
#include "Datatype.h"
#include "Op.h"
#include "Environment.h"

#endif
