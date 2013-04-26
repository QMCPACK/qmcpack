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

#ifndef _OOMPI_REQUEST_H_
#define _OOMPI_REQUEST_H_

class OOMPI_Status;
class OOMPI_Request_array;

#include <mpi.h>
#include "Status.h"
#include "Wrapper_ptr.cct"
#include "Message.h"


class OOMPI_Request
{
public:

  //
  // Constructor/Destructors
  //

  OOMPI_Request(MPI_Request a = MPI_REQUEST_NULL);
  OOMPI_Request(const OOMPI_Request &a);
  OOMPI_Request &operator=(const OOMPI_Request &a);
  OOMPI_Request &operator=(const MPI_Request &a);
  ~OOMPI_Request(void);

  //
  // Single Request functions
  //

  OOMPI_Status Wait(void);

  OOMPI_Status Test(bool &flag);
  bool Test(OOMPI_Status &s);

  void Start(void);
  void Cancel(void);
  void Free(void);

  bool Is_null(void);
  bool operator==(const OOMPI_Request &a);
  bool operator!=(const OOMPI_Request &a);

  //
  // Compatibilty
  //

  MPI_Request& Get_mpi(void);

protected:
  MPI_Request mpi_request;

private:
};


class OOMPI_Request_array
{
public:

  //
  // Constructor/Destructors
  //

  OOMPI_Request_array(int num = 1);
  OOMPI_Request_array(MPI_Request a[], int num);
  OOMPI_Request_array(const OOMPI_Request_array &a);
  OOMPI_Request_array &operator=(const OOMPI_Request_array &a);
  ~OOMPI_Request_array(void);

  //
  // Operators
  //

  OOMPI_Request& operator[](int i);

  //
  // Multiple Request_array functions
  //

  OOMPI_Status Waitany(int &index);
  int Waitany(OOMPI_Status &status);

  OOMPI_Status_array Waitall(void);
  void Waitall(OOMPI_Status_array &status);

  OOMPI_Status_array Waitsome(int &outcount, int array_of_indices[]);
  int Waitsome(OOMPI_Status_array &status, int array_of_indices[]);

  OOMPI_Status Testany(int &index, bool &flag);
  bool Testany(OOMPI_Status &status, int &index);

  OOMPI_Status_array Testall(bool &flag);
  bool Testall(OOMPI_Status_array &status);

  OOMPI_Status_array Testsome(int &outcount, int array_of_indices[]);
  int Testsome(OOMPI_Status_array &status, int array_of_indices[]);

  void Startall(void);

  void Freeall(void);

  bool operator==(const OOMPI_Request_array &a);
  bool operator!=(const OOMPI_Request_array &a);

  //
  // Compatibility
  //

  MPI_Request *Get_mpi(void);
  void Set_mpi(MPI_Request a[], int size);
  int Get_size(void);
  bool Set_size(int newsize);

protected:
  OOMPI_Request *request_array;
  MPI_Request *cache;
  int size, max_size;
  bool cache_valid, request_array_valid;

private:
  void Rebuild_request_array();
  void Rebuild_cache();

private:
};

#endif
