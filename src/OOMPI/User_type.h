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
// User_Type datatype base class
//

#ifndef _OOMPI_USER_TYPE_H_
#define _OOMPI_USER_TYPE_H_

#include <mpi.h>
#include "Message.h"


//
// This class is only used for inheritence
//

class OOMPI_User_type : public OOMPI_Message
{
public:

  inline OOMPI_User_type(OOMPI_Datatype& type, int tag)
    : OOMPI_Message(type, this, 1, tag)
  {};

  // This function exists for hysterical rasins; "top_param" parameter
  // is not necessary
  inline OOMPI_User_type(OOMPI_Datatype& type, void *top_param, int tag)
    : OOMPI_Message(type, top_param, 1, tag)
  {};

  virtual ~OOMPI_User_type(void);
};


#endif
