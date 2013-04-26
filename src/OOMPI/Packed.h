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
// Packed message class
//


#ifndef _OOMPI_PACKED_H_
#define _OOMPI_PACKED_H_

#include <mpi.h>
#include "Constants.h"
#include "Comm.h"
#include "Message.h"
#include "Datatype.h"
#include "Tag.h"


//
// Forward references
//

class OOMPI_Comm;


//
// Note: Access to MPI_Pack_size() is not provided through this object.
// It is provided through the OOMPI_Comm object, because:
//
// 1) It doesn't make sense to call MPI_Pack_size() if you already have
//    allocated a buffer to pack into/unpack from
// 2) MPI_Pack_size() takes a communicator as an argument, which is
//    hidden in the OOMPI_Comm object
//

class OOMPI_Packed : public OOMPI_Tag
{
  friend class OOMPI_Message;
  friend class OOMPI_Port;

public:

  //
  // Note that copy constructor and assignment operator are
  // both deep copies
  //

  OOMPI_Packed(int size, OOMPI_Comm &c,
               int tag = OOMPI_PACKED_TAG);
  OOMPI_Packed(void *ptr, int size, OOMPI_Comm &c,
               int tag = OOMPI_PACKED_TAG);
  OOMPI_Packed(const OOMPI_Packed &a);
  OOMPI_Packed &operator=(const OOMPI_Packed &a);
  virtual ~OOMPI_Packed();

  //
  // Access functions to get/set buffer size
  //

  int Get_size(void);
  int Set_size(int size);

  //
  // Access functions to get/set position in the buffer
  //

  int Get_position(void);
  int Set_position(int size);

  //
  // Reset the state to the beginning of the buffer
  //

  void Reset(void);

  //
  // Do the un/packing
  //

  void Start(int position = 0);

  void Pack(OOMPI_Message data);
  OOMPI_Packed &operator<<(OOMPI_Message data);
  void Pack(OOMPI_Array_message data, int count);

  void Unpack(OOMPI_Message data);
  OOMPI_Packed &operator>>(OOMPI_Message data);
  void Unpack(OOMPI_Array_message data, int count);

  void End(void);

private:

  // comm is the OOMPI_Comm that we are operating on
  // position and indicates the state of the pack
  // size is the user-specified size of the buffer
  // max_size is the actual size of the buffer (size might be < max_size)
  // buf is a pointer to the actual buffer
  // buf_created keeps track whether the buffer should be deleted
  // upon destruction or not

  MPI_Comm *comm;
  int position, size, max_size;
  char *buf;
  bool buf_created;
};

#endif
