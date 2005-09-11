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
// Status class
// 

#include <mpi.h>
#include "Status.h"
#include "Constants.h"
#include "Comm_world.h"

//
// ======================================================================
//
// OOMPI_Status
//
// ======================================================================
//

//
// MPI Constructor/ Default Constructor
//


OOMPI_Status::OOMPI_Status()
{
}


OOMPI_Status::OOMPI_Status(MPI_Status a) 
: mpi_status(a)
{
}


OOMPI_Status::OOMPI_Status(const OOMPI_Status &a)
: mpi_status(a.mpi_status)
{}


OOMPI_Status &
OOMPI_Status::operator=(const OOMPI_Status &a)
{
  if (this != &a) {
    mpi_status = a.mpi_status;
  }
  return *this;
}

OOMPI_Status &
OOMPI_Status::operator=(const MPI_Status &a)
{
  mpi_status = a;
  return *this;
}

OOMPI_Status::~OOMPI_Status(void) {}

//
// Member functions
//


// Get count of a datatype

int
OOMPI_Status::Get_count(OOMPI_Datatype datatype)
{
  int count;
  MPI_Datatype type = datatype.Get_mpi();
  MPI_Get_count(&mpi_status, type, &count);
  return count;
}

// 
// Get the element of that datatype
//

int
OOMPI_Status::Get_elements(OOMPI_Datatype datatype)
{
  int elements(0);
  MPI_Datatype type = datatype.Get_mpi();
  MPI_Get_elements(&mpi_status, type, &elements);
  return elements;
}

//
// Get source
//

int
OOMPI_Status::Get_source(void) { return mpi_status.MPI_SOURCE; }

//
// Get tag
//

int
OOMPI_Status::Get_tag(void) { return mpi_status.MPI_TAG; }

//
// Get error
//

int
OOMPI_Status::Get_error(void) { return mpi_status.MPI_ERROR; }

//
// Test to see if a status has been cancelled
//

bool
OOMPI_Status::Test_cancelled(void)
{
  int flag;
  MPI_Test_cancelled(&mpi_status, &flag);
  return (bool) flag;
}


//
// ======================================================================
//
// OOMPI_Status_array
//
// ======================================================================
//


//
// Constructors//Destructors
//

OOMPI_Status_array::OOMPI_Status_array(int num)
: size(num), max_size(num)
{
  if (size < 1) 
    size = 1;
  status_array = new OOMPI_Status [size];
  if (!status_array) 
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
  mpi_status_array = new MPI_Status [size];
  if (!mpi_status_array) 
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
  mpi_status_up_to_date = false;
}


OOMPI_Status_array::OOMPI_Status_array(MPI_Status a[], int num)
: size(0), max_size(0)
{
  Set_mpi(a, num);
}


OOMPI_Status_array::OOMPI_Status_array(const OOMPI_Status_array &a)
//: size(a.size), max_size(a.max_size)
: size(0), max_size(0)
{
#if 0
  int i;
  if (size < 1) size = 1;

  status_array = new OOMPI_Status [size];
  if (!status_array) 
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
 
  for (i = 0; i < size; i++)
    status_array[i] = a.status_array[i];
#else
  OOMPI_Status_array &src = (OOMPI_Status_array &) a;
  MPI_Status *raw_data=src.Get_mpi();
  Set_mpi(raw_data, src.Get_size());
#endif
}	 



OOMPI_Status_array &
OOMPI_Status_array::operator=(const OOMPI_Status_array &a)
{
  if (this != &a) {
    
    if (size != a.size) 
      Set_size(a.size);
    
    for (int i = 0; i < size; i++) {
      status_array[i] = a.status_array[i];
      mpi_status_array[i] = a.mpi_status_array[i];
    }

    mpi_status_up_to_date = a.mpi_status_up_to_date;
   
  }
  return *this;
}


OOMPI_Status_array::~OOMPI_Status_array(void)
{
  if (status_array)
    delete[] status_array;

  if (mpi_status_array)
    delete[] mpi_status_array;
}


//
// Operators
//

OOMPI_Status &
OOMPI_Status_array::operator[](int i)
{
  if (i < 0 || i >= size) {
    MPI_Status mpi_temp;
    mpi_temp.MPI_SOURCE = OOMPI_UNDEFINED;
    mpi_temp.MPI_TAG = OOMPI_UNDEFINED;

    static OOMPI_Status *temp =0;
    if (temp) 
      *temp = mpi_temp;
    else {
      temp = new OOMPI_Status(mpi_temp);
      if (!temp) 
	OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
    }

    //static OOMPI_Status *temp = NULL;
    //if (temp == NULL) {
    //  temp = new OOMPI_Status(mpi_temp);
    //  if (!temp) 
//	OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
 //   }
  //  else
   //   *temp = mpi_temp;

    return *temp;
  }
  return status_array[i];
}

MPI_Status *
OOMPI_Status_array::Get_mpi(void)
{
  if( !mpi_status_up_to_date ) {
    for (int i=0; i<size; i++) {
      mpi_status_array[i] = status_array[i].Get_mpi();
    }
    mpi_status_up_to_date = true;
  }

  return mpi_status_array;
}

void
OOMPI_Status_array::Set_mpi(MPI_Status a[], int num)
{
  int i;

  Set_size(num);
  
  for (i = 0; i < size; i++) {
    status_array[i] = a[i];
    mpi_status_array[i] = a[i];
  }

  mpi_status_up_to_date = true;
}	 


bool
OOMPI_Status_array::Set_size(int newsize)
{
  if (newsize == size)
    return true;

  if (newsize < 1) 
    newsize = 1;

  if (max_size < newsize) {
    OOMPI_Status *new_status_array = new OOMPI_Status [newsize];

    if (!new_status_array) {
      OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
      return false;
    }

    for (int i = 0; i < size; i++) 
      new_status_array[i] = status_array[i];

    if (size > 0) {
      delete[] status_array;
      delete[] mpi_status_array;
    }

    max_size = newsize;
    size = newsize;

    status_array = new_status_array;
    mpi_status_array = new MPI_Status [size];
    mpi_status_up_to_date = false;

    return true;
  }
  else {
    size = newsize;
    return true;
  }
}


