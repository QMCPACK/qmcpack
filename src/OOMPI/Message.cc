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
// Message datatype base class
//

#include "oompi-config.h"

#include <mpi.h>

#include "Datatype.h"
#include "Message.h"
#include "Packed.h"


//
// Local variables
//

#if OOMPI_HAVE_LONG_DOUBLE
MPI_Datatype OOMPI_Message::long_double_type;
#endif
#if OOMPI_HAVE_LONG_LONG_INT
MPI_Datatype OOMPI_Message::long_long_int_type;
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
MPI_Datatype OOMPI_Message::complex_float_type;
MPI_Datatype OOMPI_Message::complex_double_type;
MPI_Datatype OOMPI_Message::complex_long_double_type;
#endif


//
// Initialization
//
#if OOMPI_HAVE_LONG_LONG_INT || OOMPI_HAVE_LONG_DOUBLE || OOMPI_HAVE_ANSI_COMPLEX
void
OOMPI_Message::Init(void);
{
  static initialized = 0;

  if (!initialized) {
    // Take the easy way to build complex datatypes.  :)

#if OOMPI_HAVE_LONG_LONG_INT
    // JMS What to do here?
#endif
#if OOMPI_HAVE_LONG_DOUBLE
    // JMS What to do here?
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
    static OOMPI_Datatype oompi_fc;
    float_complex fc;

    oompi_fc.Contiguous(&fc, 2);
    complex_float_type= oompi_fc.Get_mpi();

    static OOMPI_Datatype oompi_dc;
    double_complex dc;

    oompi_dc.Contiguous(&dc, 2);
    complex_double_type= oompi_dc.Get_mpi();

    static OOMPI_Datatype oompi_ldc;
    long_double_complex ldc;

    oompi_ldc.Contiguous(&ldc, 2);
    complex_long_double_type= oompi_ldc.Get_mpi();
#endif

    intialized= 1;
  }
}
#endif


//
// Packed type (special)
//

OOMPI_Message::OOMPI_Message(OOMPI_Packed& i)
  : OOMPI_Tag(i.Get_tag()),
    type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_PACKED), 
    wrapped(false), top(i.buf), count(i.size)
{}

OOMPI_Message::OOMPI_Message(OOMPI_Packed& i, int tag)
  : OOMPI_Tag(i.Get_tag()),
    type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_PACKED), 
    wrapped(false), top(i.buf), count(i.size)
{
  if (tag != OOMPI_NO_TAG)
    Set_tag(tag);
}


