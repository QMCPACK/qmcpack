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


#ifndef _OOMPI_DATATYPE_H_
#define _OOMPI_DATATYPE_H_

#include <mpi.h>
#include "oompi-config.h"
#include "Wrapper_ptr.cct"
#include "Constants.h"
#include "Tag.h"
#include "Util.h"
#include "Linked_list.h"
/* #include <set> */

//
// Forward references
//

class OOMPI_Message;
class OOMPI_Array_message;
class OOMPI_Packed;
class OOMPI_Port;
class OOMPI_Datatype;
class OOMPI_Hidden;
class OOMPI_Comm_world;


//
// Class definition
//

class OOMPI_Datatype : public OOMPI_Tag
{
  friend class OOMPI_Message;
  friend class OOMPI_Array_message;
  friend class OOMPI_Comm_world;
  friend int OOMPI_INTERNAL_Free_mpi_datatype(MPI_Datatype *);

public:

  OOMPI_Datatype(MPI_Datatype a = MPI_DATATYPE_NULL,
                 int tag = OOMPI_MPI_DATATYPE_TAG);
  OOMPI_Datatype(const OOMPI_Datatype& a);
  OOMPI_Datatype& operator=(const OOMPI_Datatype& a);
  virtual ~OOMPI_Datatype(void);

  //
  // Constructors so that we can get the datatype of an unknown object
  //

  OOMPI_Datatype(OOMPI_Message& a);
  OOMPI_Datatype(OOMPI_Array_message& a);

  //
  // Constructor so that we can make static basic datatypes
  //

  OOMPI_Datatype(MPI_Datatype a, int tag, OOMPI_Hidden& hidden);

  //
  // Get the type
  //

  inline MPI_Datatype& Get_mpi(void)
  {
    return type_wrapper->Get();
  };

  //
  // MPI_DATATYPE_NULL?
  //

  inline bool Is_null(void)
  {
    return (bool) (type_wrapper->Get() == MPI_DATATYPE_NULL);
  };

  //
  // Have we been built already?
  //

  bool Built(void);

  //
  // Build simple datatypes
  //

  void Contiguous(OOMPI_Message type, int count);
  void Contiguous(OOMPI_Array_message type, int count);
  void Contiguous_type(OOMPI_Datatype type, int count);

  void Vector(int blocklength, int stride,
              OOMPI_Message type, int count);
  void Vector(int blocklength, int stride,
              OOMPI_Array_message type, int count);
  void Vector_type(int blocklength, int stride,
                   OOMPI_Datatype type, int count);

  void Hvector(int blocklength, int stride,
               OOMPI_Message type, int count);
  void Hvector(int blocklength, int stride,
               OOMPI_Array_message type, int count);
  void Hvector_type(int blocklength, int stride,
                    OOMPI_Datatype type, int count);

  void Indexed(int blocklengths[], int disps[],
               OOMPI_Message type, int count);
  void Indexed(int blocklengths[], int disps[],
               OOMPI_Array_message type, int count);
  void Indexed_type(int blocklengths[], int disps[],
                    OOMPI_Datatype type, int count);

  void Hindexed(int blocklengths[], OOMPI_Aint disps[],
                OOMPI_Message type, int count);
  void Hindexed(int blocklengths[], OOMPI_Aint disps[],
                OOMPI_Array_message type, int count);
  void Hindexed_type(int blocklengths[], OOMPI_Aint disps[],
                     OOMPI_Datatype type, int count);

  //
  // Get type info
  //

  OOMPI_Aint Extent(void);
  int Size(void);
  OOMPI_Aint Lb(void);
  OOMPI_Aint Ub(void);

  //
  // Building a structure datatype
  //

  void Struct_start(void *t, void *lb = 0);

  OOMPI_Datatype& Entry(OOMPI_Message d);
  OOMPI_Datatype& Entry(OOMPI_Array_message d, int count);
  OOMPI_Datatype& operator<<(OOMPI_Message d);

  void Struct_end(void *ub = 0);

protected:

  //
  // Internal functions and variables
  // Used for building user-defined datatypes
  //

  void Reset(void);
  MPI_Aint top;

  //
  // Reference counting wrapper
  //

  OOMPI_Wrapper_ptr<MPI_Datatype> type_wrapper;

  //
  // The set of all datatypes created.
  //

  static OOMPI_Linked_list all_datatypes;
  /* static std::set<MPI_Datatype> all_datatypes; */

  //
  // A method to free all of the datatypes in all_datatypes
  //

  static void Free_all_datatypes(void);


private:
  void inline Release(void)
  {
    OOMPI_Util a;
    a.Release_sem(SEM_DATATYPE);
  };

  //
  // Linked list of datatype info
  //

  class Type_info
  {
  public:
    Type_info(void);
    Type_info(int blocklen, MPI_Aint disp, MPI_Datatype type,
              Type_info *next);
    virtual ~Type_info();

    int MakeArrays(MPI_Aint top, int *&blocks, MPI_Aint *&disps,
                   MPI_Datatype *&types);

  private:
    int blocklen;
    MPI_Aint disp;
    MPI_Datatype type;
    Type_info *next;
  } *Type_list;

  friend class OOMPI_Datatype::Type_info;
  // To satisfy vacpp 5.0.2.0 on AIX 4.3.3 (vacpp 5.99.9999 on AIX 5.1 does
  // not complain, even without this line).  Without this line,
  // OOMPI_Datatype::Type_info::Make_arrays causes an error that
  // OOMPI_Datatype::Type_info is private.

};

#endif
