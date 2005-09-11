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

#include "Error.h"
#include "Packed.h"
#include "Comm_world.h"
#include "Constants.h"


//
// Constructors
//

OOMPI_Packed::OOMPI_Packed(int s, OOMPI_Comm &c, int t)
: OOMPI_Tag(t), comm(&c.comm_wrapper->Get()), position(0), size(s), 
  max_size(s), buf(0), buf_created(true)
{

  if (size >= 0)
    buf= new char[size];
  else
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_ARG);
}


OOMPI_Packed::OOMPI_Packed(void *ptr, int s, OOMPI_Comm &c, 
			   int t)
: OOMPI_Tag(t), comm(&c.comm_wrapper->Get()), position(0), size(s), 
  max_size(s), buf(0), buf_created(false)
{
  if (s >= 0)
    buf= (char *) ptr;
  else
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_ARG);
}


//
// Copy constructor
//

OOMPI_Packed::OOMPI_Packed(const OOMPI_Packed &a)
: OOMPI_Tag(a.Get_tag()), comm(a.comm), position(0), size(a.size),
  max_size(a.size), buf(new char[size]), buf_created(true)
{
  int i;

  if (buf != 0)
    for (i= 0; i < size; i++)
      buf[i]= a.buf[i];
  else {
    size= 0;
    max_size= 0;
    buf_created= false;
  }
}


OOMPI_Packed &
OOMPI_Packed::operator=(const OOMPI_Packed &a)
{
  int i;

  if (this != &a) {
    max_size= size= a.size;
    Set_tag(a.Get_tag());
    if (buf_created)
      delete[] buf;

    buf_created= true;
    buf= new char[size];
    if (buf != 0)
      for (i= 0; i < size; i++)
	buf[i]= a.buf[i];
    else
      max_size= size= 0;
    
    Reset();
  }

  return *this;
}


OOMPI_Packed::~OOMPI_Packed()
{
  if (buf_created)
    delete[] buf;
}


//
// Get/Set buffer size
//
int 
OOMPI_Packed::Get_size()
{
  return size;
}


int
OOMPI_Packed::Set_size(int newsize)
{
  char *newbuf;
  int i;

  // Sillyness check

  if (newsize == size)
    return newsize;

  // Check for error

  if (newsize <= 0) {
    OOMPI_ERROR.Handler(*comm, OOMPI_ERR_ARG);
    return size;
  }

  // Check to see if my allocated space can already handle it

  if (newsize <= max_size)
    size= newsize;

  // Nope, alloc more, transfer the old data, and delete the original

  else {
    newbuf= new char[newsize];
    for (i= 0; i < size; i++)
      newbuf[i]= buf[i];

    max_size= size= newsize;
    delete[] buf;
    buf= newbuf;
  }

  return size;
}


//
// Get / set position
//

int
OOMPI_Packed::Get_position()
{
  return position;
}


int
OOMPI_Packed::Set_position(int p)
{
  if (p < size && p > 0)
    position= p;

  return position;
}


//
// Reset the pack state
//

void
OOMPI_Packed::Reset()
{
  position= 0;
}


//
// Pack and unpack
//

void 
OOMPI_Packed::Start(int pos)
{
  Reset();
  position= pos;
}


void
OOMPI_Packed::Pack(OOMPI_Message data)
{
  MPI_Pack(data.Get_top(), data.Get_count(), data.Get_type(), buf, 
	   size, &position, *comm);
}


OOMPI_Packed &
OOMPI_Packed::operator<<(OOMPI_Message data)
{
  MPI_Pack(data.Get_top(), data.Get_count(), data.Get_type(), buf, 
	   size, &position, *comm);

  return *this;
}


void
OOMPI_Packed::Pack(OOMPI_Array_message data, int cnt)
{
  MPI_Pack(data.Get_top(), cnt, data.Get_type(), buf, 
	   size, &position, *comm);
}


void
OOMPI_Packed::Unpack(OOMPI_Message data)
{
  MPI_Unpack(buf, size, &position, data.Get_top(), data.Get_count(),
	     data.Get_type(), *comm);
}


OOMPI_Packed &
OOMPI_Packed::operator>>(OOMPI_Message data)
{
  MPI_Unpack(buf, size, &position, data.Get_top(), data.Get_count(),
	     data.Get_type(), *comm);

  return *this;
}


void
OOMPI_Packed::Unpack(OOMPI_Array_message data, int cnt)
{
  MPI_Unpack(buf, size, &position, data.Get_top(), cnt,
	     data.Get_type(), *comm);
}


void 
OOMPI_Packed::End(void)
{
  // Not really necessary, just for symmetry

  Reset();
}
