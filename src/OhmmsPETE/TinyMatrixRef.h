//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_TINYMATRIXREF_H
#define OHMMS_TINYMATRIXREF_H

// include files
#include "PETE/PETE.h"
#include "OhmmsPETE/TinyMatrix.h"
//////////////////////////////////////////////////////////////////////
//
// Definition of class TinyMatrixRef.
//
//////////////////////////////////////////////////////////////////////

template<class T>
class TinyMatrixRef
{
public:
  typedef T Type_t;
  TinyMatrixRef(T* datain, int d1, int d2):data(datain), D1(d1), D2(d2) { }

  template<unsigned N, unsigned M>
  TinyMatrixRef(TinyMatrix<T, N, M>& mat):
    data(mat.begin()),D1(mat.nrow()), D2(mat.ncol()) { }

  ~TinyMatrixRef()
  {
    data = NULL;
  }

  inline int nrow() const
  {
    return D1;
  }
  inline int ncol() const
  {
    return D2;
  }
  inline int byteSize() const
  {
    return D1*D2*sizeof(T);
  }

  // Get and Set Operations
  inline Type_t& operator[](unsigned int i)
  {
    return data[i];
  }
  inline Type_t operator[](unsigned int i) const
  {
    return data[i];
  }
  inline Type_t& operator()(unsigned int i)
  {
    return data[i];
  }
  inline Type_t operator()( unsigned int i) const
  {
    return data[i];
  }

  inline Type_t operator()( unsigned int i,  unsigned int j ) const
  {
    return data[i*D2+j];
  }

  inline Type_t& operator()( unsigned int i, unsigned int j )
  {
    return data[i*D2+j];
  }

private:

  TinyMatrixRef() { } // no defualt constructor is allowed
  int D1, D2;
  T* data;
};

#endif // OHMMS_TINYMATRIXREF_H

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
