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
// Group class
//

#ifndef _OOMPI_GROUP_H_
#define _OOMPI_GROUP_H_


#include <mpi.h>
#include "Comm.h"
#include "Error.h"
#include "Constants.h"


//
// Class declaraction
//

class OOMPI_Group
{
  friend class OOMPI_Comm_world;
  friend class OOMPI_Comm;
  friend class OOMPI_Intra_comm;
  friend class OOMPI_Inter_comm;

public:

  //
  // Big 4
  //

  OOMPI_Group(MPI_Group a = MPI_GROUP_NULL);
  OOMPI_Group(const OOMPI_Group &a);
  OOMPI_Group &operator=(const OOMPI_Group &a);
  inline ~OOMPI_Group() {};

  //
  // Group accessors
  //

  // returns the size of the group

  int Size(void);

  // returns the rank of the calling process

  int Rank(void);

  // access MPI_Group

  inline MPI_Group &Get_mpi(void)
  {
    return group_wrapper->Get();
  };

  void Translate_ranks(int n, int ranks1[], OOMPI_Group g2, int ranks2[]);
  int *Translate_ranks(int n, int ranks1[], OOMPI_Group g2);

  OOMPI_Group Incl(int n, int ranks[]);
  OOMPI_Group Excl(int n, int ranks[]);
  OOMPI_Group Range_incl(int n, int ranges[][3]);
  OOMPI_Group Range_excl(int n, int ranges[][3]);

  //
  // Group comparison operations
  //

  OOMPI_Compare Compare(OOMPI_Group g2);

  // compare g1 and g2: returns TRUE if MPI_EQUAL or MPI_SIMILAR
  friend bool operator== (OOMPI_Group g1, OOMPI_Group g2);

  // compare g1 and g2: returns TRUE if MPI_UNEQUAL
  friend bool operator!= (OOMPI_Group g1, OOMPI_Group g2);

  //
  // Group set operations
  //

  // take the union of g1 and g2: returns a new group
  friend OOMPI_Group operator| (OOMPI_Group g1, OOMPI_Group g2);

  // take the intersection of g1 and g2: returns a new group
  friend OOMPI_Group operator& (OOMPI_Group g1, OOMPI_Group g2);

  // take the difference of g1 and g2: - returns a new group
  friend OOMPI_Group operator- (OOMPI_Group g1, OOMPI_Group g2);

  // Inquiries

  bool Is_null(void);
  bool Is_empty(void);

protected:
  OOMPI_Wrapper_ptr<MPI_Group> group_wrapper;

private:
};


//
// Global instances
//

extern const OOMPI_Group OOMPI_GROUP_NULL;
extern const OOMPI_Group OOMPI_GROUP_EMPTY;

#endif
