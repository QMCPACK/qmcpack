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
// Message datatype class variable
//

//
// Forward reference
//

class OOMPI_Datatype;

#include <mpi.h>
#include "Comm.h"
#include "Comm_world.h"
#include "Message.h"
#include "Datatype.h"
#include "Error.h"
#include "Hidden.h"

//
// Make all_datatypes
//

// JJW We should replace this with STL when we do the STL stuff:
// std::set<MPI_Datatype> OOMPI_Datatype::all_datatypes;
OOMPI_Linked_list OOMPI_Datatype::all_datatypes;

//
// Wrapper delete function
// This has to have the strange name because it is external (required
// to make it a friend function).
//
int 
OOMPI_INTERNAL_Free_mpi_datatype(MPI_Datatype *ptr)
{

  MPI_Datatype temp_datatype = *ptr;

  if (ptr != 0 && *ptr != MPI_DATATYPE_NULL) {
    // We only need to free if we haven't already finalized.
    if (!OOMPI_COMM_WORLD.Finalized())
      MPI_Type_free(ptr);
  }

  OOMPI_Datatype::all_datatypes.erase(temp_datatype);

  return MPI_SUCCESS;
}



//
// Constructor (doubles as MPI constructor)
//
OOMPI_Datatype::OOMPI_Datatype(MPI_Datatype a, int tag)
  : OOMPI_Tag(tag), type_wrapper(a, 0), Type_list(0)
{
  // We don't have to clean these up
}


//
// Copy constructor
//
OOMPI_Datatype::OOMPI_Datatype(const OOMPI_Datatype &a) 
: OOMPI_Tag(a.my_tag), type_wrapper(a.type_wrapper), Type_list(0)
{
}


//
// Assignment operator
//
OOMPI_Datatype &
OOMPI_Datatype::operator=(const OOMPI_Datatype &a)
{
  Type_list = 0;

  if (this != &a) {
    type_wrapper = a.type_wrapper;
    my_tag = a.my_tag;
  }

  return *this;
}


//
// Virtual destructor
//
OOMPI_Datatype::~OOMPI_Datatype()
{
}


//
// Special constructors
//
OOMPI_Datatype::OOMPI_Datatype(OOMPI_Message &m)
  : OOMPI_Tag(m.Get_tag()), Type_list(0)
{
  if (m.wrapped)
    type_wrapper = m.type_wrapper;
  else
    type_wrapper->Set(m.native_type, 0);
}


OOMPI_Datatype::OOMPI_Datatype(OOMPI_Array_message &m)
  : OOMPI_Tag(m.Get_tag()), Type_list(0)
{
  // It can only be a basic datatype -- it will not be wrapped!
  type_wrapper->Set(m.Get_type(), 0);
}


//
// Constructor so that we can make static basic datatypes
//

OOMPI_Datatype::OOMPI_Datatype(MPI_Datatype a, int tag, OOMPI_Hidden& hidden)
  : OOMPI_Tag(tag), Type_list(0)
{
  type_wrapper->Set(a, 0);
  hidden.do_nothing(); // For picky compilers
}


//
// Built
//
bool
OOMPI_Datatype::Built()
{
  OOMPI_Util a;
  a.Get_sem(SEM_DATATYPE);
  if (!Is_null()) {
    a.Release_sem(SEM_DATATYPE);
    return true;
  }
  return false;
}


//
// MPI_Type_contiguous
//
void
OOMPI_Datatype::Contiguous(OOMPI_Message type, int cnt)
{
  Contiguous_type(type.Get_type(), cnt);
}


void
OOMPI_Datatype::Contiguous(OOMPI_Array_message type, int cnt)
{
  Contiguous_type(type.Get_type(), cnt);
}


void
OOMPI_Datatype::Contiguous_type(OOMPI_Datatype type, int cnt)
{
  MPI_Datatype *out= new MPI_Datatype;

  Reset();
  if (MPI_Type_contiguous(cnt, type.type_wrapper->Get(), out) != MPI_SUCCESS) {
    Release();
    return;
  }

  if (MPI_Type_commit(out) != MPI_SUCCESS) {
    Release();
    return;
  }
  
  type_wrapper->Set(*out, OOMPI_INTERNAL_Free_mpi_datatype);
  all_datatypes.insert(Get_mpi());
  delete out;
  Release();
}


void
OOMPI_Datatype::Vector(int blocklength, int stride, 
		       OOMPI_Message type, int cnt)
{
  Vector_type(blocklength, stride, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Vector(int blocklength, int stride, 
		       OOMPI_Array_message type, int cnt)
{
  Vector_type(blocklength, stride, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Vector_type(int blocklength, int stride, 
			    OOMPI_Datatype type, int cnt)
{
  MPI_Datatype *out= new MPI_Datatype;

  Reset();
  if (MPI_Type_vector(cnt, blocklength, stride, type.type_wrapper->Get(),
		      out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }
  if (MPI_Type_commit(out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }

  type_wrapper->Set(*out, OOMPI_INTERNAL_Free_mpi_datatype);
  all_datatypes.insert(Get_mpi());
  delete out;
  Release();
}


void
OOMPI_Datatype::Hvector(int blocklength, int stride, 
			OOMPI_Message type, int cnt)
{
  Hvector_type(blocklength, stride, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Hvector(int blocklength, int stride, 
			OOMPI_Array_message type, int cnt)
{
  Hvector_type(blocklength, stride, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Hvector_type(int blocklength, int stride, 
			     OOMPI_Datatype type, int cnt)
{
  MPI_Datatype *out= new MPI_Datatype;

  Reset();
  if (MPI_Type_hvector(cnt, blocklength, stride, 
		       type.type_wrapper->Get(), out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }
  if (MPI_Type_commit(out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }

  type_wrapper->Set(*out, OOMPI_INTERNAL_Free_mpi_datatype);
  all_datatypes.insert(Get_mpi());
  delete out;
  Release();
}


void
OOMPI_Datatype::Indexed(int blocklengths[], int disps[], 
			OOMPI_Message type, int cnt)
{
  Indexed_type(blocklengths, disps, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Indexed(int blocklengths[], int disps[], 
			OOMPI_Array_message type, int cnt)
{
  Indexed_type(blocklengths, disps, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Indexed_type(int blocklengths[], int disps[], 
			     OOMPI_Datatype type, int cnt)
{
  MPI_Datatype *out= new MPI_Datatype;
  
  Reset();
  if (MPI_Type_indexed(cnt, blocklengths, disps, 
		       type.type_wrapper->Get(), out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }
  if (MPI_Type_commit(out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }

  type_wrapper->Set(*out, OOMPI_INTERNAL_Free_mpi_datatype);
  all_datatypes.insert(Get_mpi());
  delete out;
  Release();
}


void
OOMPI_Datatype::Hindexed(int blocklengths[], OOMPI_Aint disps[], 
			 OOMPI_Message type, int cnt)
{
  Hindexed_type(blocklengths, disps, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Hindexed(int blocklengths[], OOMPI_Aint disps[], 
			 OOMPI_Array_message type, int cnt)
{
  Hindexed_type(blocklengths, disps, type.Get_type(), cnt);
}


void
OOMPI_Datatype::Hindexed_type(int blocklengths[], OOMPI_Aint disps[], 
			      OOMPI_Datatype type, int cnt)
{
  MPI_Datatype *out= new MPI_Datatype;

  Reset();
  if (MPI_Type_hindexed(cnt, blocklengths, (MPI_Aint *) disps, 
			type.type_wrapper->Get(), out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }
  if (MPI_Type_commit(out) != MPI_SUCCESS) {
    Release();
    delete out;
    return;
  }

  type_wrapper->Set(*out, OOMPI_INTERNAL_Free_mpi_datatype);
  all_datatypes.insert(Get_mpi());
  delete out;
  Release();
}


OOMPI_Aint
OOMPI_Datatype::Extent()
{
  MPI_Aint extent;

  if (type_wrapper->Get() == 0) {
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_ARG);
    return OOMPI_UNDEFINED;
  }

  if (MPI_Type_extent(type_wrapper->Get(), &extent) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return (OOMPI_Aint) extent;
}


int
OOMPI_Datatype::Size()
{
  int size;

  if (type_wrapper->Get() == 0) {
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_ARG);
    return OOMPI_UNDEFINED;
  }

  if (MPI_Type_size(type_wrapper->Get(), &size) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return size;
}


OOMPI_Aint 
OOMPI_Datatype::Lb()
{
  MPI_Aint lb;

  if (type_wrapper->Get() == 0) {
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_ARG);
    return OOMPI_UNDEFINED;
  }

  if (MPI_Type_lb(type_wrapper->Get(), &lb) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return (OOMPI_Aint) lb;
}


OOMPI_Aint
OOMPI_Datatype::Ub()
{
  MPI_Aint ub;

  if (type_wrapper->Get() == 0) {
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_ARG);
    return OOMPI_UNDEFINED;
  }

  if (MPI_Type_lb(type_wrapper->Get(), &ub) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return (OOMPI_Aint) ub;
}


//
// Building a structure datatype
//

void
OOMPI_Datatype::Struct_start(void *t, void *lb)
{
  top = (MPI_Aint) t;
  Reset();

  if (lb != 0) {
    MPI_Aint disp;

    MPI_Address(lb, &disp);
    Type_list = new Type_info(1, disp, MPI_LB, Type_list);
  }
}


OOMPI_Datatype &
OOMPI_Datatype::Entry(OOMPI_Message d)
{
  MPI_Aint disp;

  MPI_Address(d.Get_top(), &disp);
  Type_list = new Type_info(d.Get_count(), disp, d.Get_type(), Type_list);

  return *this;
}


OOMPI_Datatype &
OOMPI_Datatype::Entry(OOMPI_Array_message d, int this_count)
{
  MPI_Aint disp;

  MPI_Address(d.Get_top(), &disp);
  Type_list = new Type_info(this_count, disp, d.Get_type(), Type_list);

  return *this;
}


OOMPI_Datatype &
OOMPI_Datatype::operator<<(OOMPI_Message d)
{
  MPI_Aint disp;

  MPI_Address(d.Get_top(), &disp);
  Type_list = new Type_info(d.Get_count(), disp, d.Get_type(), Type_list);

  return *this;
}


void
OOMPI_Datatype::Struct_end(void *ub)
{
  int *blocks = (int *) 0, cnt = 0;
  MPI_Aint *disps = (MPI_Aint *) 0, disp;
  MPI_Datatype *types = (MPI_Datatype *) 0, *newtype;

  // Sillyness check

  if (Type_list == 0) {
    Release();
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_TYPE);
    return;
  }

  // Did the user supply an upper bound?

  if (ub != 0) {
    MPI_Address(ub, &disp);
    Type_list = new Type_info(1, disp, MPI_UB, Type_list);
  }

  // Build the datatype

  cnt = Type_list->MakeArrays(top, blocks, disps, types);
  delete Type_list;

  newtype = new MPI_Datatype;
  if (MPI_Type_struct(cnt, blocks, disps, types, newtype) != MPI_SUCCESS) {
    Release();
    return;
  }
  if (MPI_Type_commit(newtype) != MPI_SUCCESS)  {
    Release();
    return;
  }

  type_wrapper->Set(*newtype, OOMPI_INTERNAL_Free_mpi_datatype);
  all_datatypes.insert(Get_mpi());
  delete newtype;
  delete[] blocks;
  delete[] disps;
  delete[] types;
  Release();
}


/////////////////////////////////////////////////


//
// Reset the state of this instance so that we can start building
// a new datatype
//
void 
OOMPI_Datatype::Reset()
{
  Set_tag(OOMPI_MESSAGE_TAG);

  if (Type_list != 0)
    delete Type_list;
  Type_list = 0;

  type_wrapper->Unset();
}


void OOMPI_Datatype::Free_all_datatypes(void)
{
  all_datatypes.deleteAll();
  all_datatypes.clear();
  #if 0
  std::set<MPI_Datatype>::const_iterator i;
  for(i = all_datatypes.begin(); i != all_datatypes.end(); i++) {
    MPI_Datatype temp_datatype;
    temp_datatype = *i;
    MPI_Type_free(&temp_datatype);
  }
  #endif /* Enable for STL */
}

/////////////////////////////////////////////////


//
// Type_info
// Default constructor
//
OOMPI_Datatype::Type_info::Type_info()
{
  blocklen = 0;
  disp = 0;
  type = MPI_DATATYPE_NULL;
  next = 0;
}


//
// Type_info
// Constructor
//
OOMPI_Datatype::Type_info::Type_info(int b, MPI_Aint d, MPI_Datatype t,
				     Type_info *n)
{
  blocklen = b;
  disp = d;
  type = t;
  next = n;
}


//
// Type_info
// Destructor
//
OOMPI_Datatype::Type_info::~Type_info()
{
  if (next != 0)
    delete next; // Hey, that's cool!  ;-)
}


//
// Type_info
// MakeArrays
//
int
OOMPI_Datatype::Type_info::MakeArrays(MPI_Aint top, int *&blocks, 
				      MPI_Aint *&disps, 
				      MPI_Datatype *&types)
{
  int cnt, i;
  OOMPI_Datatype::Type_info *cur;

  // Count how many we need
  // Also get the first offset.  IMPORTANT: The *first* offset
  // is the *last* member of the list!

  for (cnt= 0, cur= this; cur != 0; cnt++, cur = cur->next)
    continue;

  // Sillyness check

  if (cnt == 0) {
    blocks = 0;
    disps = 0;
    types = 0;

    return 0;
  }

  // Alloc out the space

  blocks = new int[cnt];
  disps = new MPI_Aint[cnt];
  types = new MPI_Datatype[cnt];

  // Fill them in

  for (i = cnt - 1, cur = this; cur != 0; i--, cur = cur->next) {
    blocks[i] = cur->blocklen;
    disps[i] = cur->disp - top;
    types[i] = cur->type;
  }
  
  return cnt;
}

