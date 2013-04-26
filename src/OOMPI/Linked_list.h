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
// Linked list of MPI datatypes (for cleanup at finalize)
//

#ifndef _OOMPI_LINKED_LIST_H_
#define _OOMPI_LINKED_LIST_H_

class OOMPI_Linked_list_node;

class OOMPI_Linked_list
{
  // This class defines a linked list of MPI datatypes, and a method to
  // delete all of them.  This should be replaced by an STL set when we
  // get STL enabled in OOMPI.  The list is unsorted, and all list
  // algorithms are the standard ones (new items go at the beginning of
  // the list).

public:
  // Create on OOMPI_Linked_list object.
  OOMPI_Linked_list();

  // Delete the linked list.
  ~OOMPI_Linked_list();

  // Insert a datatype into the list.
  void insert(MPI_Datatype data);

  // Remove a datatype from the list.
  void erase(MPI_Datatype data);

  // MPI_Type_free all types in the list.
  void deleteAll();

  // Delete the contents of the list (without MPI_Type_free'ing them).
  void clear();

private:
  OOMPI_Linked_list_node *head;
};

struct OOMPI_Linked_list_node
{
  // This class contains one node of the linked list.  This is for internal
  // use of the OOMPI_Linked_list class only.

  MPI_Datatype data;
  OOMPI_Linked_list_node *next;
};

#endif /* _OOMPI_LINKED_LIST_H_ */
