//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.   
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc. 
//////////////////////////////////////////////////////////////////////////////////////


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
