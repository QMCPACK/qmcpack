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


#ifndef _OOMPI_STATUS_H_
#define _OOMPI_STATUS_H_


//
// Forward reference
//

#include <mpi.h>
#include "Wrapper_ptr.cct"
#include "Datatype.h"
#include "Message.h"

class OOMPI_Status_array;

//
// OOMPI_Status
//

class OOMPI_Status
{
public:

  //
  // Constructors//Destructors
  //

  OOMPI_Status(void);
  OOMPI_Status(MPI_Status a);
  OOMPI_Status(const OOMPI_Status &a);
  OOMPI_Status &operator=(const OOMPI_Status &a);
  OOMPI_Status &operator=(const MPI_Status &a);
  ~OOMPI_Status(void);

  //
  // Access functions
  //

  int Get_count(OOMPI_Datatype datatype);
  int Get_elements(OOMPI_Datatype datatype);
  int Get_tag(void);
  int Get_source(void);
  int Get_error(void);
  bool Test_cancelled(void);

  //
  // Compatibility
  //

  inline MPI_Status &Get_mpi(void)
  {
    return mpi_status;
  }

protected:
  MPI_Status mpi_status;

private:

};


//
// OOMPI_Status_array
//


class OOMPI_Status_array
{
public:

  //
  // Constructors//Destructors
  //

  OOMPI_Status_array(int num = 1);
  OOMPI_Status_array(MPI_Status a[], int num);
  OOMPI_Status_array(const OOMPI_Status_array &a);
  OOMPI_Status_array& operator=(const OOMPI_Status_array &a);
  ~OOMPI_Status_array(void);

  //
  // Operators
  //

  OOMPI_Status & operator[](int i);

  //
  // Compatibility
  //

  MPI_Status *Get_mpi(void);
  void Set_mpi(MPI_Status a[], int num);
  inline int Get_size(void)
  {
    return size;
  };
  bool Set_size(int newsize);

protected:
  OOMPI_Status *status_array;
  MPI_Status *mpi_status_array;
  bool mpi_status_up_to_date;
  int size;
  int max_size;

private:
};

#endif








