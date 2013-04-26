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
// Message datatype base class
//

#ifndef _OOMPI_MESSAGE_H_
#define _OOMPI_MESSAGE_H_

#include "oompi-config.h"

#include <mpi.h>
#if OOMPI_HAVE_ANSI_COMPLEX
#include <complex>
#endif

#include "Wrapper_ptr.cct"
#include "Datatype.h"
#include "Constants.h"
#include "Tag.h"


class OOMPI_Message : public OOMPI_Tag
{
  friend class OOMPI_Datatype;

public:

  //
  // Direct constructor
  //

  inline OOMPI_Message(const OOMPI_Datatype& type,
                       void* d, int cnt = 1)
    : OOMPI_Tag(type.Get_tag()),
      type_wrapper(type.type_wrapper),
      native_type(MPI_DATATYPE_NULL),
      wrapped(true), top(d), count(cnt)
  {}

  inline OOMPI_Message(const OOMPI_Datatype& type, void* d,
                       int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(type.type_wrapper),
      native_type(MPI_DATATYPE_NULL),
      wrapped(true), top(d), count(cnt)
  {}

  //
  // Copy constructor
  //

  inline OOMPI_Message(const OOMPI_Message& m)
    : OOMPI_Tag(m.Get_tag()),
      type_wrapper(m.type_wrapper),
      native_type(m.native_type),
      wrapped(m.wrapped), top(m.top), count(m.count)
  {};

  //
  // Assignment operator
  //

  inline OOMPI_Message& operator=(OOMPI_Message& m)
  {
    if (&m != this)
    {
      type_wrapper = m.type_wrapper;
      native_type = m.native_type;
      wrapped = m.wrapped;
      top = m.top;
      count = m.count;
      Set_tag(m.Get_tag());
    }
    return *this;
  };

  //
  // Destructor
  //

  inline ~OOMPI_Message() {};

  //
  // Initialization (build supplemental datatypes)
  //

#if OOMPI_HAVE_LONG_LONG_INT || OOMPI_HAVE_LONG_DOUBLE || OOMPI_HAVE_ANSI_COMPLEX
  void Init(void);
#endif

  //
  // Implicit promotions from base types
  //

  inline OOMPI_Message(char& i)
    : OOMPI_Tag(OOMPI_CHAR_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_CHAR),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(short& i)
    : OOMPI_Tag(OOMPI_SHORT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_SHORT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(int& i)
    : OOMPI_Tag(OOMPI_INT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_INT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long& i)
    : OOMPI_Tag(OOMPI_LONG_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_LONG),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned char& i)
    : OOMPI_Tag(OOMPI_UNSIGNED_CHAR_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_CHAR),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned short& i)
    : OOMPI_Tag(OOMPI_UNSIGNED_SHORT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_SHORT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned& i)
    : OOMPI_Tag(OOMPI_UNSIGNED_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned long& i)
    : OOMPI_Tag(OOMPI_UNSIGNED_LONG_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_LONG),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(float& i)
    : OOMPI_Tag(OOMPI_FLOAT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_FLOAT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(double& i)
    : OOMPI_Tag(OOMPI_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_DOUBLE),
      wrapped(false), top(&i), count(1)
  {}

#if OOMPI_HAVE_LONG_LONG_INT
  inline OOMPI_Message(long long int& i)
    : OOMPI_Tag(OOMPI_LONG_LONG_INT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), type(oompi_long_long_int_type),
      wrapped(false), top(&i), count(1)
  {}
#endif

#if OOMPI_HAVE_LONG_DOUBLE
  inline OOMPI_Message(long double& i)
    : OOMPI_Tag(OOMPI_LONG_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), type(long_double_type),
      wrapped(false), top(&i), count(1)
  {}
#endif

#if OOMPI_HAVE_ANSI_COMPLEX
  inline OOMPI_Message(float_complex& i)
    : OOMPI_Tag(OOMPI_COMPLEX_FLOAT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), type(complex_float_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(double_complex& i)
    : OOMPI_Tag(OOMPI_COMPLEX_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), type(complex_double_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long_double_complex& i)
    : OOMPI_Tag(OOMPI_COMPLEX_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), type(complex_long_double_type),
      wrapped(false), top(&i), count(1)
  {}
#endif


  //
  // Explicit instantiations
  //

  inline OOMPI_Message(char& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_CHAR),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(char* i, int cnt)
    : OOMPI_Tag(OOMPI_CHAR_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_CHAR),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(char* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_CHAR),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(short& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_SHORT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(short* i, int cnt)
    : OOMPI_Tag(OOMPI_SHORT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_SHORT),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(short* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_SHORT),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(int& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_INT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(int* i, int cnt)
    : OOMPI_Tag(OOMPI_INT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_INT),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(int* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_INT),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(long& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_LONG),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long* i, int cnt)
    : OOMPI_Tag(OOMPI_LONG_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_LONG),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(long* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_LONG),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(unsigned char& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_CHAR),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned char* i, int cnt)
    : OOMPI_Tag(OOMPI_UNSIGNED_CHAR_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_CHAR),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(unsigned char* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_CHAR),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(unsigned short& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_SHORT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned short* i, int cnt)
    : OOMPI_Tag(OOMPI_UNSIGNED_SHORT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_SHORT),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(unsigned short* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_SHORT),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(unsigned& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned* i, int cnt)
    : OOMPI_Tag(OOMPI_UNSIGNED_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(unsigned* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(unsigned long& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_LONG),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(unsigned long* i, int cnt)
    : OOMPI_Tag(OOMPI_UNSIGNED_LONG_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_LONG),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(unsigned long* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_UNSIGNED_LONG),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(float& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_FLOAT),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(float* i, int cnt)
    : OOMPI_Tag(OOMPI_FLOAT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_FLOAT),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(float* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_FLOAT),
      wrapped(false), top(i), count(cnt)
  {}

  inline OOMPI_Message(double& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_DOUBLE),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(double* i, int cnt)
    : OOMPI_Tag(OOMPI_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_DOUBLE),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(double* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(MPI_DOUBLE),
      wrapped(false), top(i), count(cnt)
  {}

#if OOMPI_HAVE_LONG_LONG_INT
  inline OOMPI_Message(long long int& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(long_long_int_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long long int* i, int cnt)
    : OOMPI_Tag(OOMPI_LONG_LONG_INT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(long_long_int_type),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(long long int* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(long_long_int_type),
      wrapped(false), top(i), count(cnt)
  {}
#endif

#if OOMPI_HAVE_LONG_DOUBLE
  inline OOMPI_Message(long double& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(long_double_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long double* i, int cnt)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(long_double_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long double* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(long_double_type),
      wrapped(false), top(i), count(cnt)
  {}
#endif

#if OOMPI_HAVE_ANSI_COMPLEX
  inline OOMPI_Message(float_complex* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(complex_float_type),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(float_complex& i)
    : OOMPI_Tag(OOMPI_COMPLEX_FLOAT_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(complex_float_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(float_complex& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(complex_float_type),
      wrapped(false), top(&i), count(1)
  {}

  inline OOMPI_Message(double_complex* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(complex_double_type),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(double_complex& i)
    : OOMPI_Tag(OOMPI_COMPLEX_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(complex_double_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(double_complex& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0), native_type(complex_double_type),
      wrapped(false), top(&i), count(1)
  {}

  inline OOMPI_Message(long_double_complex* i, int cnt, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0),
      native_type(complex_long_double_type),
      wrapped(false), top(i), count(cnt)
  {}
  inline OOMPI_Message(long_double_complex& i)
    : OOMPI_Tag(OOMPI_COMPLEX_LONG_DOUBLE_TAG),
      type_wrapper(MPI_DATATYPE_NULL, 0),
      native_type(complex_long_double_type),
      wrapped(false), top(&i), count(1)
  {}
  inline OOMPI_Message(long_double_complex& i, int tag)
    : OOMPI_Tag(tag),
      type_wrapper(MPI_DATATYPE_NULL, 0),
      native_type(complex_long_double_type),
      wrapped(false), top(&i), count(1)
  {}
#endif

  //
  // Packed type (special)
  // Not inlined because of #include file dependencies
  //

  OOMPI_Message(OOMPI_Packed& i);
  OOMPI_Message(OOMPI_Packed& i, int tag);


  //
  // Access functions
  //

  inline MPI_Datatype Get_type(void)
  {
    return (wrapped) ? type_wrapper->Get() : native_type;
  };
  inline void* Get_top(void)
  {
    return top;
  };
  inline int Get_count(void)
  {
    return count;
  };

protected:

  //
  // Message information
  //

  OOMPI_Wrapper_ptr<MPI_Datatype> type_wrapper;
  MPI_Datatype native_type;
  bool wrapped;

  void* top;
  int count;

private:
};


//
// Same class, but for arrays
//

class OOMPI_Array_message : public OOMPI_Tag
{
public:
  //
  // Promotions from base types
  //

  inline OOMPI_Array_message(char i[])
    : OOMPI_Tag(OOMPI_CHAR_TAG), type(MPI_CHAR), top(i)
  {}
  inline OOMPI_Array_message(short i[])
    : OOMPI_Tag(OOMPI_SHORT_TAG), type(MPI_SHORT), top(i)
  {}
  inline OOMPI_Array_message(int i[])
    : OOMPI_Tag(OOMPI_INT_TAG), type(MPI_INT), top(i)
  {}
  inline OOMPI_Array_message(long i[])
    : OOMPI_Tag(OOMPI_LONG_TAG), type(MPI_LONG), top(i)
  {}
  inline OOMPI_Array_message(unsigned char i[])
    : OOMPI_Tag(OOMPI_UNSIGNED_CHAR_TAG), type(MPI_UNSIGNED_CHAR), top(i)
  {}
  inline OOMPI_Array_message(unsigned short i[])
    : OOMPI_Tag(OOMPI_UNSIGNED_SHORT_TAG), type(MPI_UNSIGNED_SHORT), top(i)
  {}
  inline OOMPI_Array_message(unsigned i[])
    : OOMPI_Tag(OOMPI_UNSIGNED_TAG), type(MPI_UNSIGNED), top(i)
  {}
  inline OOMPI_Array_message(unsigned long i[])
    : OOMPI_Tag(OOMPI_UNSIGNED_LONG_TAG), type(MPI_UNSIGNED_LONG), top(i)
  {}
  inline OOMPI_Array_message(float i[])
    : OOMPI_Tag(OOMPI_FLOAT_TAG), type(MPI_FLOAT), top(i)
  {}
  inline OOMPI_Array_message(double i[])
    : OOMPI_Tag(OOMPI_DOUBLE_TAG), type(MPI_DOUBLE), top(i)
  {}

#if OOMPI_HAVE_LONG_LONG_INT
  inline OOMPI_Array_message(long long int i[])
    : OOMPI_Tag(OOMPI_LONG_LONG_INT_TAG), type(long_long_int_type), top(i)
  {}
#endif

#if OOMPI_HAVE_LONG_DOUBLE
  inline OOMPI_Array_message(long double i[])
    : OOMPI_Tag(OOMPI_LONG_DOUBLE_TAG), type(long_double_type), top(i)
  {}
#endif

#if OOMPI_HAVE_ANSI_COMPLEX
  inline OOMPI_Array_message(float_complex i[])
    : OOMPI_Tag(OOMPI_COMPLEX_FLOAT_TAG), type(complex_float_type), top(i)
  {}
  inline OOMPI_Array_message(double_complex i[])
    : OOMPI_Tag(OOMPI_COMPLEX_DOUBLE_TAG), type(complex_double_type), top(i)
  {}
  inline OOMPI_Array_message(long_double_complex i[])
    : OOMPI_Tag(OOMPI_COMPLEX_LONG_DOUBLE_TAG),
      type(complex_long_double_type), top(i)
  {}
#endif

  //
  // Access functions
  //

  inline MPI_Datatype Get_type(void)
  {
    return type;
  };
  inline void* Get_top(void)
  {
    return top;
  };

protected:

  //
  // Message information
  //

  MPI_Datatype type;
  void* top;

  //
  // Class variables
  //

#if OOMPI_HAVE_LONG_DOUBLE
  static MPI_Datatype long_double_type;
#endif
#if OOMPI_HAVE_LONG_LONG_INT
  static MPI_Datatype long_long_int_type;
#endif
#if OOMPI_HAVE_ANSI_COMPLEX
  static MPI_Datatype complex_float_type;
  static MPI_Datatype complex_double_type;
  static MPI_Datatype complex_long_double_type;
#endif

private:
};

#endif
