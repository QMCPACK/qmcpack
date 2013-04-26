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
// Cartesian communicators
//

#ifndef _OOMPI_CART_COMM_H_
#define _OOMPI_CART_COMM_H_

#include "oompi-config.h"
#include "Comm.h"
#include "Intra_comm.h"
#include <stdarg.h>


class OOMPI_Cart_comm : public OOMPI_Intra_comm
{
  friend class OOMPI_Intra_comm;

public:

  // Big 4

  // default constructor.  sets comm to MPI_COMM_NULL
  inline OOMPI_Cart_comm(MPI_Comm mpi_comm = MPI_COMM_NULL);

  // Copy constructor.  Does a shallow copy
  OOMPI_Cart_comm(const OOMPI_Cart_comm &a);

  // Assignment operator.  Does a shallow copy
  OOMPI_Cart_comm& operator=(const OOMPI_Cart_comm &);

  // Creates a Cart_comm from an Intra_comm
  // Does an MPI_Cart_create
  OOMPI_Cart_comm(OOMPI_Intra_comm& intra_comm, int ndims,
                  int dims[], bool periods[], bool reorder = false);

  // If the dimensions are all either periodic, or not periodic, a
  // single bool can be used.
  OOMPI_Cart_comm(OOMPI_Intra_comm& intra_comm, int ndims,
                  int dims[], bool periods, bool reorder = false);


  // Destructor
  ~OOMPI_Cart_comm();

  // ---------- Communicator management

  // MPI_Comm_dup
  OOMPI_Cart_comm Dup(void);

  // ---------- Process Topologies

  // MPI_Cart_sub
  OOMPI_Cart_comm Sub(bool remain_dims[]);
  OOMPI_Cart_comm Sub(bool dim0, ...);

  // MPI_Cartdim_get
  int Dim_get(void);

  // MPI_Get
  // if maxdims is not given then it defaults to the number of dims
  void Get(int maxdims, int dims[], bool period[] = 0, int coords[] = 0);
  void Get(int dims[], bool period[] = 0, int coords[] = 0);

  // MPI_Cart_rank
  int Rank(int coords[]);
  int Rank(int, ...);
  int Rank(void);

  // MPI_Cart_rank (returns a port)
  OOMPI_Port operator() (int coords[]);
  OOMPI_Port operator() (int, ...);

  // MPI_Cart_coords
  // if coords is not supplied, the function allocated memory.
  // if maxdims is not specified the it used all the dims.
  // if rank is not specified it uses the rank of the calling process.
  void Coords(int rank, int maxdims, int coords[]);
  void Coords(int rank, int coords[]);
  void Coords(int coords[]);
  int *Coords(int rank, int maxdims);
  int *Coords(int rank);
  int *Coords(void);

  //  MPI_Cart_shift
  int Shift(int direction, int disp, int &src);
  int Shift(int direction, int disp);

  // return flag
  virtual bool Test_inter(void);

protected:
private:
  // Used for various array assemblies and conversions (bool -> int)
  int *int_array;

  // Internal constructors (default + needs_to_be_freed)
  void do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed);
  // This exists since constructors can't be called directly.
  inline OOMPI_Cart_comm(MPI_Comm mpi_comm, bool needs_to_be_freed);
};

// MPI and default constructor
inline OOMPI_Cart_comm::OOMPI_Cart_comm(MPI_Comm mpi_comm)
  : OOMPI_Intra_comm(MPI_COMM_NULL), int_array(0)
{
  do_full_init(mpi_comm, false);
}

// Full constructor
inline OOMPI_Cart_comm::OOMPI_Cart_comm(MPI_Comm mpi_comm, bool needs_to_be_freed)
  : OOMPI_Intra_comm(MPI_COMM_NULL), int_array(0)
{
  do_full_init(mpi_comm, needs_to_be_freed);
}

#endif
