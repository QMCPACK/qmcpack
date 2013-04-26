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
// Communicator base class
//

#ifndef _OOMPI_HIDDEN_H_
#define _OOMPI_HIDDEN_H_

#include "Op.h"

//
// This is a "hidden" class so that we can have some private
// constructors.
//

class OOMPI_Hidden
{
public:
  static inline OOMPI_Op Create_op(MPI_Op op)
  {
    OOMPI_Hidden bogus;
    return OOMPI_Op(op, bogus);
  }

  static inline OOMPI_Datatype Create_datatype(MPI_Datatype type, int tag)
  {
    OOMPI_Hidden bogus;
    return OOMPI_Datatype(type, tag, bogus);
  }

  static inline void do_nothing(void) {};
protected:
  int i;
private:
  // We never want an instance of this object

  inline OOMPI_Hidden() {};
};

#endif
