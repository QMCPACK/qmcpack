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
// Request class
//

#include "oompi-config.h"
#include "Request.h"
#include "Status.h"
#include "Comm_world.h"


//
// ======================================================================
//
// OOMPI_Request
//
// ======================================================================
//

OOMPI_Request::OOMPI_Request(MPI_Request a) 
: mpi_request(a)
{
}
 
 
OOMPI_Request::OOMPI_Request(const OOMPI_Request &a)
: mpi_request(a.mpi_request)
{
}

 
OOMPI_Request &
OOMPI_Request::operator=(const OOMPI_Request &a)
{
  if (this != &a)
    mpi_request = a.mpi_request;
  return *this;
}

 
OOMPI_Request &
OOMPI_Request::operator=(const MPI_Request &a)
{
  mpi_request = a;
  return *this;
}

 
OOMPI_Request::~OOMPI_Request(void) 
{
}


//
// Wait access function
//

OOMPI_Status 
OOMPI_Request::Wait(void)
{
  MPI_Status mpi_status;
  MPI_Wait(&Get_mpi(), &mpi_status);
  return OOMPI_Status(mpi_status);
}


//
// Test access function
//

OOMPI_Status
OOMPI_Request::Test(bool &flag)
{
  MPI_Status mpi_status;
#if OOMPI_BOOL_NE_INT
  int intflag;
  MPI_Test(&Get_mpi(), &intflag, &mpi_status);
  flag = (bool) intflag;
#else
  MPI_Test(&Get_mpi(), (int *) &flag, &mpi_status);
#endif
  return OOMPI_Status(mpi_status);
}


bool
OOMPI_Request::Test(OOMPI_Status &s)
{
  int flag;
  MPI_Test(&(Get_mpi()), &flag, &(s.Get_mpi() ));
  return (bool) flag;
}

//
// Start access function
//

void
OOMPI_Request::Start(void)
{
  MPI_Start(&Get_mpi());
}

//
// Cancel access function
//

void
OOMPI_Request::Cancel(void)
{
  MPI_Cancel(&Get_mpi());
}

//
// Free this request
//

void
OOMPI_Request::Free(void)
{
  // JJW Is_null should be inline
#if !OOMPI_IBM21014
  // SP MPI does not seem to like this free.  It gives an "internal error"
  // when -euidevelop is on, and acts strangely when it is off.
  if (!Is_null())
    MPI_Request_free(&Get_mpi());
#endif
}

//
// Is MPI_REQUEST_NULL
//

bool
OOMPI_Request::Is_null(void)
{
  return (bool) (mpi_request == MPI_REQUEST_NULL);
}

bool
OOMPI_Request::operator==(const OOMPI_Request &a)
{
  return (bool) (mpi_request == a.mpi_request);
}

bool
OOMPI_Request::operator!=(const OOMPI_Request &a)
{
  return (bool) (mpi_request != a.mpi_request);
}

//
// Get_mpi
//

MPI_Request &
OOMPI_Request::Get_mpi()
{
  return mpi_request;
}


//
// ======================================================================
//
// OOMPI_Request_array
//
// ======================================================================
//
 
 
//
// Constructors//Destructors
//
 

OOMPI_Request_array::OOMPI_Request_array(int num)
: size(num), max_size(num), cache_valid(true), request_array_valid(true)
{
  int i;
  
  if (size < 1) {
    size = 1;
    max_size = 1;
  }

  request_array = new OOMPI_Request[size];
  cache = new MPI_Request [size];

  for (i=0; i<size; i++)
    request_array[i] = cache[i] = MPI_REQUEST_NULL;
}
 

OOMPI_Request_array::OOMPI_Request_array(MPI_Request a[], int num)
: size(0), max_size(0), request_array_valid(true)
{
  Set_mpi(a, num);
}


OOMPI_Request_array::OOMPI_Request_array(const OOMPI_Request_array &a)
//: size(a.size), cache_valid(true), request_array_valid(true), max_size(a.size)
: size(0), max_size(0), request_array_valid(true)
{
#if 0
  int i;

  if (size < 1) {
    size = 1;
    max_size = 1;
  }

  request_array = new OOMPI_Request [size];
  cache = new MPI_Request [size];

  if (a.request_array_valid) {
    for (i = 0; i < size; i++) {
      request_array[i] = a.request_array[i];
      cache[i] = a.request_array[i].Get_mpi();
    }
  }
  else {
    for (i = 0; i < size; i++) {
      request_array[i] = a.cache[i];
      cache[i] = a.cache[i];
    }
  }
#else
  OOMPI_Request_array &src = (OOMPI_Request_array &) a;
  Set_mpi(src.Get_mpi(), src.Get_size());
#endif
}        
 
 
 
OOMPI_Request_array &
OOMPI_Request_array::operator=(const OOMPI_Request_array &a)
{
  if (this != &a) {

    if (size != a.size) 
      Set_size(a.size);
	
    if (a.request_array_valid == true) {
      for (int i=0; i<size; i++) {
	request_array[i] = a.request_array[i];
	cache[i] = a.request_array[i].Get_mpi();
      }
    }
    else {
      for (int i=0; i<size; i++) {
	request_array[i] = a.cache[i];
	cache[i] = a.cache[i];
      }
    }
    cache_valid = true;
    request_array_valid = true;
  }
  return *this;
}

OOMPI_Request_array::~OOMPI_Request_array(void)
{
  if (request_array)
    delete[] request_array;

  if (cache)
    delete[] cache;
}


//
// Operators
//
 

OOMPI_Request &
OOMPI_Request_array::operator[](int i)
{
  if (i < 0 || i >= size) {
    //AIX error, complains NULL
    static OOMPI_Request *temp = 0;
    if(temp)
      *temp = MPI_REQUEST_NULL;
    else
      temp = new OOMPI_Request(MPI_REQUEST_NULL);
    //static OOMPI_Request *temp = NULL;
    //if (temp == NULL) 
    //  temp = new OOMPI_Request(MPI_REQUEST_NULL);
    //else
    //  *temp = MPI_REQUEST_NULL;
 
    return *temp;
  }
  if (!request_array_valid) 
    Rebuild_request_array();
  cache_valid = false;
  return request_array[i];
}


OOMPI_Status
OOMPI_Request_array::Waitany(int &index)
{
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Status mpi_status;
  MPI_Waitany(size, cache, &index, &mpi_status);
  request_array_valid = false;

  return OOMPI_Status(mpi_status);
}

int
OOMPI_Request_array::Waitany(OOMPI_Status &a)
{
  int index;
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Waitany(size, cache, &index, &(a.Get_mpi()));
  request_array_valid = false;

  return index;
}

OOMPI_Status_array
OOMPI_Request_array::Waitall(void)
{
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Status *mpi_status = new MPI_Status [size];
  MPI_Waitall(size, cache, mpi_status);
  request_array_valid = false;

  OOMPI_Status_array S(mpi_status, size);
  delete[] mpi_status;

  return S;
}

void 
OOMPI_Request_array::Waitall(OOMPI_Status_array &a)
{
  if (!cache_valid) 
    Rebuild_cache();

  if (a.Get_size() != size)
    a.Set_size(size);
  MPI_Status *mpi_status = a.Get_mpi();
  MPI_Waitall(size, cache, mpi_status);
  a.Set_mpi(mpi_status, size);

  request_array_valid = false;
}

OOMPI_Status_array
OOMPI_Request_array::Waitsome(int &outcount, int array_of_indices[])
{
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Status *mpi_status = new MPI_Status [size];
  MPI_Waitsome(size, cache, &outcount, array_of_indices, mpi_status);

  request_array_valid = false;

  OOMPI_Status_array S(mpi_status, size);
  delete[] mpi_status;

  return S;
}

int 
OOMPI_Request_array::Waitsome(OOMPI_Status_array &a, int array_of_indices[])
{
  int outcount;
  if (!cache_valid) 
    Rebuild_cache();

  if (a.Get_size() != size)
    a.Set_size(size);
  MPI_Status *mpi_status = a.Get_mpi();
  MPI_Waitsome(size, cache, &outcount, array_of_indices, mpi_status);
  a.Set_mpi(mpi_status, size);

  request_array_valid = false;
  return outcount;
}

OOMPI_Status
OOMPI_Request_array::Testany(int &index, bool &flag)
{
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Status mpi_status;
#if OOMPI_BOOL_NE_INT
  int intflag;
  MPI_Testany(size, cache, &index, &intflag, &mpi_status);
  flag = (bool) intflag;
#else
  MPI_Testany(size, cache, &index, (int *) &flag, &mpi_status);
#endif

  request_array_valid = false;
  return OOMPI_Status(mpi_status);
}

bool
OOMPI_Request_array::Testany(OOMPI_Status &a, int &index)
{
  if (!cache_valid) 
    Rebuild_cache();

  int flag;
  MPI_Testany(size, cache, &index, &flag, &(a.Get_mpi()));

  request_array_valid = false;

  return (bool) flag;
}

OOMPI_Status_array
OOMPI_Request_array::Testall(bool &flag)
{
  if (!cache_valid)
    Rebuild_cache();
  MPI_Status *mpi_status = new MPI_Status [size];

#if OOMPI_BOOL_NE_INT
  int intflag;
  MPI_Testall(size, cache, &intflag, mpi_status);
  flag = (bool) intflag;
#else
  MPI_Testall(size, cache, (int *) &flag, mpi_status);
#endif

  request_array_valid = false;
  OOMPI_Status_array S(mpi_status, size);
  delete[] mpi_status;

  return S;
}

bool
OOMPI_Request_array::Testall(OOMPI_Status_array &a)
{
  if (!cache_valid)
    Rebuild_cache();

  if (a.Get_size() != size)
    a.Set_size(size);
  MPI_Status *mpi_status = a.Get_mpi();

  int intflag;
  MPI_Testall(size, cache, &intflag, mpi_status);
  a.Set_mpi(mpi_status, size);

  request_array_valid = false;

  return ((!intflag) ? false : true);
}

OOMPI_Status_array 
OOMPI_Request_array::Testsome(int &outcount, int array_of_indices[])
{
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Status *mpi_status = new MPI_Status [size];
  MPI_Testsome(size, cache, &outcount, array_of_indices, mpi_status);

  request_array_valid = false;
  OOMPI_Status_array S(mpi_status, size);
  delete[] mpi_status;

  return S;
}

int 
OOMPI_Request_array::Testsome(OOMPI_Status_array &a, int array_of_indices[])
{
  int outcount;
  if (!cache_valid) 
    Rebuild_cache();

  if (a.Get_size() != size)
    a.Set_size(size);
  MPI_Status *mpi_status = a.Get_mpi();
  MPI_Testsome(size, cache, &outcount, array_of_indices, mpi_status);
  a.Set_mpi(mpi_status, size);

  request_array_valid = false;
  return outcount;
}

void
OOMPI_Request_array::Startall()
{
  if (!cache_valid) 
    Rebuild_cache();
  MPI_Startall(size, cache);
  request_array_valid = false;
}

void
OOMPI_Request_array::Freeall()
{
  if (!request_array_valid)
    Rebuild_request_array();
  for(int i=0; i<size; i++)
    request_array[i].Free();
}


MPI_Request *
OOMPI_Request_array::Get_mpi()
{
  if (!cache_valid) 
    Rebuild_cache();
  return cache;
}

void
OOMPI_Request_array::Set_mpi(MPI_Request a[], int num)
{
  int i;

  Set_size(num);

  for (i = 0; i < size; i++)
    request_array[i] = cache[i] = a[i];

  cache_valid = true;
  request_array_valid = true;
}        
 
int 
OOMPI_Request_array::Get_size(void)
{
  return size;
}

bool
OOMPI_Request_array::Set_size(int newsize)
{
  if (!request_array_valid) 
    Rebuild_request_array();
  
  if (newsize == size)
    return true;

  if (newsize < 1) 
    newsize = 1;

  if (max_size < newsize) {
    OOMPI_Request *new_request_array = new OOMPI_Request [newsize];
    MPI_Request *new_cache = new MPI_Request [newsize];

    if (!new_request_array || !new_cache)
      OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
    
    for (int i = 0; i < size; i++) {
      new_request_array[i] = request_array[i];
      new_cache[i] = request_array[i].Get_mpi();
    }

    if (size > 0) {
      delete[] request_array;
      delete[] cache;
    }

    max_size = newsize;
    size = newsize;

    request_array = new_request_array;
    cache = new_cache;
    return true;
  }
  
  else {
    size = newsize;
    return true;
  }
}

bool
OOMPI_Request_array::operator==(const OOMPI_Request_array &a)
{
  if (size != a.size)
    return false;
  for (int i=0; i<size; i++)
    if (request_array[i] != a.request_array[i])
      return false;
  
  return true;
}

bool
OOMPI_Request_array::operator!=(const OOMPI_Request_array &a)
{
  return (bool) !(*this == a);
}

void
OOMPI_Request_array::Rebuild_cache()
{
  for (int i = 0; i < size; i++) 
    cache[i] = request_array[i].Get_mpi();
  cache_valid = true;
}

void
OOMPI_Request_array::Rebuild_request_array()
{
  for (int i = 0; i < size; i++) 
    if (request_array[i].Get_mpi() != cache[i])
      request_array[i] = cache[i];
  request_array_valid = true;
}
