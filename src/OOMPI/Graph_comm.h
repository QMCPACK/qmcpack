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
// Intracommunicators
//


#ifndef _OOMPI_GRAPH_COMM_H_
#define _OOMPI_GRAPH_COMM_H_

#include "Group.h"
#include "Comm.h"
#include "Intra_comm.h"


class OOMPI_Graph_comm : public OOMPI_Intra_comm
{
  friend class OOMPI_Intra_comm;

public:

  // Big 4

  // default constructor.  sets comm to MPI_COMM_NULL
  inline OOMPI_Graph_comm(MPI_Comm mpi_comm = MPI_COMM_NULL);

  // Shallow copies
  OOMPI_Graph_comm(const OOMPI_Graph_comm& a);
  OOMPI_Graph_comm& operator=(const OOMPI_Graph_comm& a);

  // Does an MPI_Graph_create
  OOMPI_Graph_comm(OOMPI_Intra_comm& intra_comm, int nnodes,
                   int index[], int edges[], bool reorder = false);

  // Destructor
  ~OOMPI_Graph_comm(void);

  // ---------- Communicator management

  // Does an MPI_Comm_dup
  OOMPI_Graph_comm Dup(void);

  // ---------- Process Topologies

  // MPI_Graphdims_get
  void Dims_get(int& nnodes, int& nedges);
  int Num_nodes(void);
  int Num_edges(void);

  // MPI_Graph_get
  // if 0 is used for maxindex, or maxedges then Num_edges() or
  // Num_nodes will be used.
  void Get(int maxindex, int maxedges, int index[], int edges[]);
  void Get(int index[], int edges[]);
  void Get_index(int index[]);
  void Get_edges(int edges[]);
  int* Get_index(void);
  int* Get_edges(void);


  // MPI_Graph_Neighbors
  // Where neighbors is not provided, the function will allocate the
  // memory. Where maxneighbors is not provided, MPI_Graph_neighbors_count
  // will be used.  Where rank is not provided, the rank of the current
  // process is used.

  int* Neighbors(int maxneighbors, int rank, int neighbors[]);
  int* Neighbors(int rank, int neighbors[]);
  int* Neighbors(int neighbors[]);
  int* Neighbors(int maxneighbors, int rank);
  int* Neighbors(int rank);
  int* Neighbors(void);

  // MPI_Graph_neighbors_count
  // where rank is not provided, it returns the number of neighbors of
  // the current process.

  int Neighbors_count(int rank);
  int Neighbors_count(void);

  // virtual function from OOMPI_Comm
  virtual bool Test_inter(void);

protected:
private:

  // Standard constructors (default + needs_to_be_freed)
  inline OOMPI_Graph_comm(MPI_Comm mpi_comm, bool needs_to_be_freed);
  void do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed);
  // This exists since constructors can't be called directly.

};

// MPI and default constructor
inline OOMPI_Graph_comm::OOMPI_Graph_comm(MPI_Comm mpi_comm)
  : OOMPI_Intra_comm(MPI_COMM_NULL)
{
  do_full_init(mpi_comm, false);
}

// Full constructor
inline OOMPI_Graph_comm::OOMPI_Graph_comm(MPI_Comm mpi_comm, bool needs_to_be_freed)
  : OOMPI_Intra_comm(MPI_COMM_NULL)
{
  do_full_init(mpi_comm, needs_to_be_freed);
}

#endif

