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

#include <mpi.h>
#include "Linked_list.h"

// Create on OOMPI_Linked_list object.
OOMPI_Linked_list::OOMPI_Linked_list() 
{
  head = 0;
}

// Delete the linked list.
OOMPI_Linked_list::~OOMPI_Linked_list() 
{
  clear();
}

// Insert a datatype into the list.
void OOMPI_Linked_list::insert(MPI_Datatype data) 
{
  OOMPI_Linked_list_node *old_head = head;
  head = new OOMPI_Linked_list_node;
  head->data = data;
  head->next = old_head;
}

// Remove a datatype from the list.
void OOMPI_Linked_list::erase(MPI_Datatype data) 
{
  // Delete all elements of the list equal to the given data value.
  // Since the list is unsorted, all of it must be searched.  The
  // beginning of the list is a special case, since its pointer is
  // not stored in an OOMPI_Linked_list_node.

  OOMPI_Linked_list_node *old_head;
  while (head != 0 && head->data == data) {
    old_head = head;
    head = head->next;
    delete old_head;
  }

  OOMPI_Linked_list_node *current;
  OOMPI_Linked_list_node *old_next;

  for (current = head; current != 0 && current->next != 0; 
       current = current->next)
    while (current->next->data == data) {
      old_next = current->next;
      current->next = current->next->next;
      delete old_next;
    }
}

// MPI_Type_free all types in the list.
void OOMPI_Linked_list::deleteAll() 
{
  OOMPI_Linked_list_node *current;

  for (current = head; current != 0; current = current->next)
    MPI_Type_free(&current->data);
}

// Delete the contents of the list (without MPI_Type_free'ing them).
void OOMPI_Linked_list::clear() 
{
  OOMPI_Linked_list_node *old_head;

  while (head != 0) {
    old_head = head;
    head = head->next;
    delete old_head;
  }
}

