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
// Cart_communicators
// 

#include <mpi.h>
#include <stdarg.h>
#include "Comm.h"
#include "Cart_comm.h"
#include "Inter_comm.h"
#include "Error.h"


void OOMPI_Cart_comm::do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed)
{
  if (mpi_comm == MPI_COMM_NULL)
    return;

  int status(0);
  MPI_Topo_test(mpi_comm, &status);

  if (status == MPI_CART) {
    MPI_Constructor(mpi_comm, needs_to_be_freed);
    if (mpi_comm != MPI_COMM_NULL) {
      int ndims(0);
      MPI_Cartdim_get(mpi_comm, &ndims);
      int_array = new int[ndims];
    }
  }
  else 
    OOMPI_ERROR.Handler(mpi_comm, OOMPI_ERR_COMM);
}


OOMPI_Cart_comm::OOMPI_Cart_comm(const OOMPI_Cart_comm &a)
: OOMPI_Intra_comm(a), int_array(0)
{
  if (comm_wrapper->Get() != MPI_COMM_NULL) {
    int ndims(0);
    MPI_Cartdim_get(comm_wrapper->Get(), &ndims);
    int_array = new int[ndims];
  }
}


OOMPI_Cart_comm &
OOMPI_Cart_comm::operator=(const OOMPI_Cart_comm &a) 
{
  (OOMPI_Intra_comm &) *this = (OOMPI_Intra_comm &) a; 

  if (comm_wrapper->Get() != MPI_COMM_NULL) {
    int ndims(0);
    MPI_Cartdim_get(comm_wrapper->Get(), &ndims);
    if (int_array != 0)
      delete[] int_array;
    int_array = new int[ndims];
  }
  else
    int_array = 0;

  return *this;
}


// Does an MPI_Cart_create
OOMPI_Cart_comm::OOMPI_Cart_comm(OOMPI_Intra_comm &intra_comm, int ndims, 
				 int dims[], bool periods[], bool reorder)
: OOMPI_Intra_comm(MPI_COMM_NULL), int_array(0)
{
  MPI_Comm mpi_comm(0);
  int *int_ptr;
  int_array = new int[ndims];

#if OOMPI_BOOL_NE_INT
  int i;
  int_ptr = int_array;
  for (i = 0; i < ndims; i++)
    int_array[i] = (int) periods[i];
#else
  int_ptr = (int *) periods;
#endif
  if (MPI_Cart_create(intra_comm.comm_wrapper->Get(), ndims, 
		      dims, int_ptr, reorder, &mpi_comm) 
      != MPI_SUCCESS)
    mpi_comm = MPI_COMM_NULL;
  
  MPI_Constructor(mpi_comm, true);
}

OOMPI_Cart_comm::OOMPI_Cart_comm(OOMPI_Intra_comm &intra_comm, int ndims, 
				 int dims[], bool periods, bool reorder)
: OOMPI_Intra_comm(MPI_COMM_NULL), int_array(0)
{
  MPI_Comm mpi_comm(0);
  int_array = new int[ndims];
  for (int i = 0; i < ndims; i++)
    int_array[i] = (int) periods;

  if (MPI_Cart_create(intra_comm.comm_wrapper->Get(), ndims, 
		      dims, int_array, reorder, &mpi_comm) 
      != MPI_SUCCESS) {
    mpi_comm = MPI_COMM_NULL;
  }

  MPI_Constructor(mpi_comm, true);
}


OOMPI_Cart_comm::~OOMPI_Cart_comm()
{
  if (int_array != 0)
    delete[] int_array;
}


// ---------- Communicator management

// Does an MPI_Comm_dup
OOMPI_Cart_comm
OOMPI_Cart_comm::Dup(void)
{
  MPI_Comm mpi_comm;
  if (MPI_Comm_dup(Get_mpi(), &mpi_comm) != MPI_SUCCESS)
    return OOMPI_Cart_comm(MPI_COMM_NULL);

  return OOMPI_Cart_comm(mpi_comm, true);
}


// ---------- Process Topologies


// returns newcomm
OOMPI_Cart_comm
OOMPI_Cart_comm::Sub(bool remain_dims[])
{
  int *int_ptr;
  MPI_Comm mpi_comm(0);

#if OOMPI_BOOL_NE_INT
  int i, ndims = Dim_get();
  int_ptr = int_array;
  for (i = 0; i < ndims; i++)
    int_array[i] = (int) remain_dims[i];
#else
  int_ptr = (int *) remain_dims;
#endif

  if (MPI_Cart_sub(comm_wrapper->Get(), int_ptr, &mpi_comm) 
      != MPI_SUCCESS)
    return OOMPI_Cart_comm(MPI_COMM_NULL);
  
  return OOMPI_Cart_comm(mpi_comm, true);
}

OOMPI_Cart_comm
OOMPI_Cart_comm::Sub(bool dim0, ...)
{
  va_list ap;
  int dims = Dim_get();
  MPI_Comm mpi_comm(0);
 
  int_array[0] = (int) dim0;
  va_start(ap, dim0);
  for (int i = 1; i < dims; i++)
    int_array[i] = (int) va_arg(ap, int);
  va_end(ap);
  
  if (MPI_Cart_sub(comm_wrapper->Get(), int_array, &mpi_comm) 
      != MPI_SUCCESS) {
    return OOMPI_Cart_comm(MPI_COMM_NULL);
  }
  
  return OOMPI_Cart_comm(mpi_comm, true);
}

// MPI_Cartdim_get
int 
OOMPI_Cart_comm::Dim_get()
{
  int ndims(0);
  MPI_Cartdim_get(comm_wrapper->Get(), &ndims);
  return ndims;
}


// MPI_Cart_get
void
OOMPI_Cart_comm::Get(int maxdims, int dims[], bool period[], int coords[])
{
  bool free_period(false), free_coords(false);
  int *period_ptr=0, *coords_ptr=0;

  if (period == 0) {
    free_period = true;
    period_ptr = new int[maxdims];
  }
  else {
#if OOMPI_BOOL_NE_INT
    if (maxdims <= Dim_get())
      period_ptr = int_array;
    else {
      // This is an error, but let the MPI implementation handle it
      // (the user cannot request more dimensions than are in this comm).
      // At least they won't get an arbitrary memory error because of OOMPI!
      period_ptr = new int[maxdims];
      free_period = true;
    }
#else
    period_ptr = (int *) &(period[0]);
#endif
  }

  if (coords == 0) {
    free_coords = true;
    coords_ptr = new int [maxdims];
  }
  else 
    coords_ptr = &(coords[0]);

  MPI_Cart_get(comm_wrapper->Get(), maxdims, dims, period_ptr, 
	       coords_ptr);

#if OOMPI_BOOL_NE_INT
  if (!free_period) {
    int i;
    for (i = 0; i < maxdims; i++)
      period[i] = (period_ptr[i] == 1 ? true : false);
  }
#endif

  if (free_period)
    delete[] period_ptr;
  if (free_coords)
    delete[] coords_ptr;
}

// MPI_Cart_get
void
OOMPI_Cart_comm::Get(int dims[], bool period[], int coords[])
{
  Get(Dim_get(), dims, period, coords);
}


// MPI_Cart_rank
int 
OOMPI_Cart_comm::Rank(int coords[])
{
  int rank(0);
  MPI_Cart_rank(comm_wrapper->Get(), coords, &rank);
  return rank;
}

// MPI_Cart_rank with variable length arguments
int
OOMPI_Cart_comm::Rank(int dim0, ...)
{
  va_list ap;
  int dims = Dim_get();
  MPI_Comm comm = comm_wrapper->Get();

  int_array[0] = dim0;
  va_start(ap, dim0);
  for (int i = 1; i < dims; i++)
    int_array[i] = va_arg(ap, int);
  va_end(ap);
  
  int rank(0);
  if (MPI_Cart_rank(comm, int_array, &rank) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return rank;
}

int
OOMPI_Cart_comm::Rank()
{
  return OOMPI_Intra_comm::Rank();
}

// MPI_Cart_rank (returns a Port)
OOMPI_Port
OOMPI_Cart_comm::operator() (int coords[])
{
  int rank(0);
  MPI_Cart_rank(comm_wrapper->Get(), coords, &rank);
  return (*this)[rank];
}

// MPI_Cart_rank with variable length arguments (returns a Port)
OOMPI_Port
OOMPI_Cart_comm::operator() (int dim0, ...)
{
  va_list ap;
  int dims = Dim_get();
  MPI_Comm comm = comm_wrapper->Get();

  int_array[0] = dim0;
  va_start(ap, dim0);
  for (int i = 1; i < dims; i++)
    int_array[i] = va_arg(ap, int);
  va_end(ap);
  
  int rank(0);
  if (MPI_Cart_rank(comm, int_array, &rank) != MPI_SUCCESS)
    return OOMPI_Port();
  
  return (*this)[rank];
}


// MPI_Cart_coords
void
OOMPI_Cart_comm::Coords(int rank, int maxdims, int coords[])
{
  MPI_Cart_coords(comm_wrapper->Get(), rank, maxdims, coords);
  return;
}

// MPI_Cart_coords
void
OOMPI_Cart_comm::Coords(int rank, int coords[])
{
  MPI_Cart_coords(comm_wrapper->Get(), rank, Dim_get(), coords);
  return;
}

// MPI_Cart_coords
void
OOMPI_Cart_comm::Coords(int coords[])
{
  MPI_Cart_coords(comm_wrapper->Get(), OOMPI_Comm::Rank(), Dim_get(), coords);
  return;
}

// MPI_Cart_coords
int *
OOMPI_Cart_comm::Coords(int rank, int maxdims)
{
  int *coords = new int [maxdims];
  if (MPI_Cart_coords(comm_wrapper->Get(), rank, maxdims, coords) 
      != MPI_SUCCESS) {
    delete[] coords;
    return (int *) 0;
  }

  return coords;
}

// MPI_Cart_coords
int *
OOMPI_Cart_comm::Coords(int rank)
{
  int *coords = new int [Dim_get()];
  if (MPI_Cart_coords(comm_wrapper->Get(), rank, Dim_get(), coords) 
      != MPI_SUCCESS) {
    delete[] coords;
    return (int *) 0;
  }

  return coords;
}

// MPI_Cart_coords
int *
OOMPI_Cart_comm::Coords()
{
  int *coords = new int [Dim_get()];
  if (MPI_Cart_coords(comm_wrapper->Get(), OOMPI_Comm::Rank(), 
		      Dim_get(), coords) 
      != MPI_SUCCESS) {
    delete[] coords;
    return (int *) 0;
  }

  return coords;
}


// MPI_Cart_shift
int
OOMPI_Cart_comm::Shift(int direction, int disp, int &rank_source)
{
  int rank_dest(0);
  if (MPI_Cart_shift(comm_wrapper->Get(), direction, disp, 
		  &rank_source, &rank_dest) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;
  return rank_dest;
}

// MPI_Cart_shift
int
OOMPI_Cart_comm::Shift(int direction, int disp)
{
  int rank_dest(0), rank_source(0);
  if (MPI_Cart_shift(comm_wrapper->Get(), direction, disp, 
		  &rank_source, &rank_dest) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;
  return rank_dest;
}

// ----------- Type - virtual function from comm
// returns flag  
bool
OOMPI_Cart_comm::Test_inter()
{
  return OOMPI_Intra_comm::Test_inter();
}

