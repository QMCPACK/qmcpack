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

#include "oompi-config.h"

#include <iostream>

#include "Group.h"
#include "Comm_world.h"


const OOMPI_Group OOMPI_GROUP_EMPTY(MPI_GROUP_EMPTY);
const OOMPI_Group OOMPI_GROUP_NULL(MPI_GROUP_NULL);


//
// Wrapper delete function
//

static int 
Free_mpi_group(MPI_Group *ptr)
{
  if (ptr != 0 && *ptr != MPI_GROUP_NULL && *ptr != MPI_GROUP_EMPTY ) {
    if (OOMPI_COMM_WORLD.Finalized()) {
      std::cerr << "Attempt to free group after finalize (ignored, probably resulting in memory leak)." 
		<< std::endl;
    } else {
      MPI_Group_free(ptr);
    }
  }

  return MPI_SUCCESS;
}


//
// Default constructor
// Disgusting cast for picky compilers
//

OOMPI_Group::OOMPI_Group(MPI_Group a)
  : group_wrapper(a, ((a == MPI_GROUP_EMPTY || a == MPI_GROUP_NULL) ? 
		      (int (*)(MPI_Group* arg)) 0 : Free_mpi_group))
{
}


//
// Copy constructor
//

OOMPI_Group::OOMPI_Group(const OOMPI_Group &a)
  : group_wrapper(a.group_wrapper)
{
}


//
// Assignment operator
//

OOMPI_Group &
OOMPI_Group::operator=(const OOMPI_Group &a)
{
  group_wrapper = a.group_wrapper;

  return (*this);
}


//
// ---------- Group accessors ----------
//

//
// Get the size of the Group
//

int 
OOMPI_Group::Size(void)
{
  int size;
  
  if (MPI_Group_size(group_wrapper->Get(), &size) != MPI_SUCCESS)
    return -1;
  
  return size;
}


//
// Get the rank of the calling process
//

int
OOMPI_Group::Rank()
{
  int rank;

  if (MPI_Group_rank(group_wrapper->Get(), &rank) != MPI_SUCCESS)
    return -1;
  
  return rank;
}


//
// Determine the relative numbering of the same processes on two 
// different groups
//

void
OOMPI_Group::Translate_ranks(int n, int ranks1[], OOMPI_Group g2, 
			     int ranks2[])
{
  MPI_Group_translate_ranks(group_wrapper->Get(), n, ranks1, 
			    g2.group_wrapper->Get(), ranks2);
}


//
// Translate_ranks
// This version allocates memory for the ranks
//

int *
OOMPI_Group::Translate_ranks(int n, int ranks1[], OOMPI_Group g2)
{
  int *ranks2 = new int [n];
  if (MPI_Group_translate_ranks(group_wrapper->Get(), n, ranks1, 
				g2.group_wrapper->Get(), 
				ranks2) != MPI_SUCCESS) {
    delete [] ranks2;
    return (int *) 0;
  }
  return ranks2;
}


//
// Compares two groups and returns MPI_IDENT, MPI_SIMILAR, or MPI_NOT_EQUAL
//

OOMPI_Compare
OOMPI_Group::Compare(OOMPI_Group g2)
{
  int result;
  if (MPI_Group_compare(group_wrapper->Get(), 
			g2.group_wrapper->Get(), &result) != MPI_SUCCESS)
    return OOMPI_UNEQUAL;

  return (OOMPI_Compare) result;
}


//
// Creates a new group of processes consisting of ranks from 0 to n-1
//

OOMPI_Group
OOMPI_Group::Incl(int n, int ranks[])
{
  MPI_Group mpi_group;

  if (MPI_Group_incl(group_wrapper->Get(), n, ranks, &mpi_group) != 
      MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}


//
// Creates a new group of processes by deleting ranks from 0 to n-1
//

OOMPI_Group
OOMPI_Group::Excl(int n, int ranks[])
{
  MPI_Group mpi_group;

  if (MPI_Group_excl(group_wrapper->Get(), n, ranks, &mpi_group) != 
      MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}


//
// Creates a new group of processes by including a range
//

OOMPI_Group 
OOMPI_Group::Range_incl(int n, int ranges[][3])
{
  MPI_Group mpi_group;

  if (MPI_Group_range_incl(group_wrapper->Get(), n, ranges, &mpi_group) != 
      MPI_SUCCESS) 
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}

//
// Creates a new group of processes by deleting a range
//

OOMPI_Group
OOMPI_Group::Range_excl(int n, int ranges[][3])
{
  MPI_Group mpi_group;

  if (MPI_Group_range_excl(group_wrapper->Get(), n, ranges, &mpi_group) != 
      MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}


//
// ---------- Group comparison operations ----------
//

//
// compare g1 and g2: returns TRUE if NOT NOT_EQUAL
//

bool
operator== (OOMPI_Group g1, OOMPI_Group g2)
{
  return (bool) !(g1 != g2);
}


//
// compare g1 and g2: returns TRUE if MPI_UNEQUAL
//

bool
operator!= (OOMPI_Group g1, OOMPI_Group g2)
{
  int result;
  
  if (MPI_Group_compare(g1.group_wrapper->Get(), g2.group_wrapper->Get(), 
			&result) != MPI_SUCCESS)
    return true;
  return (bool) (result == MPI_UNEQUAL);
}


//
// take the union of g1 and g2: returns a new group
//

OOMPI_Group
operator| (OOMPI_Group g1, OOMPI_Group g2)
{
  MPI_Group mpi_group;
  
  if (MPI_Group_union(g1.group_wrapper->Get(), g2.group_wrapper->Get(), 
		      &mpi_group) != MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);
  
  return OOMPI_Group(mpi_group);
}


//
// take the intersection of g1 and g2: returns a new group
//

OOMPI_Group
operator& (OOMPI_Group g1, OOMPI_Group g2)
{
  MPI_Group mpi_group;
  
  if (MPI_Group_intersection(g1.group_wrapper->Get(), g2.group_wrapper->Get(),
			     &mpi_group) != MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}


//
// take the difference of g1 and g2: returns a new group 
//

OOMPI_Group
operator-(OOMPI_Group g1, OOMPI_Group g2)
{
  MPI_Group mpi_group;

  if (MPI_Group_difference(g1.group_wrapper->Get(), g2.group_wrapper->Get(), 
			   &mpi_group) != MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}

//
// Inqueries
//

bool
OOMPI_Group::Is_null(void)
{
  return (bool) (group_wrapper->Get() == MPI_GROUP_NULL);
}


bool
OOMPI_Group::Is_empty(void)
{
  return (bool) (group_wrapper->Get() == MPI_GROUP_EMPTY);
}

  
