//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_TINYMATRIX_H
#define OHMMS_TINYMATRIX_H

// include files
#include "PETE/PETE.h"
#include "OhmmsPETE/OhmmsTinyMeta.h"

//////////////////////////////////////////////////////////////////////
//
// Definition of class TinyMatrix.
//
//////////////////////////////////////////////////////////////////////

template<class T, unsigned D1, unsigned D2>
class TinyMatrix
{
public:

  typedef T Type_t;
  enum { Nrow = D1, Ncol = D2};

  // D1,D2efault Constructor initializes to zero.
  TinyMatrix()
  {
    OTAssign<TinyMatrix<T,D1,D2>, T, OpAssign>::apply(*this,T(0), OpAssign());
  }

  // A noninitializing ctor.
  class DontInitialize {};
  TinyMatrix(DontInitialize) {}

  // Copy Constructor
  TinyMatrix(const TinyMatrix<T,D1,D2> &rhs)
  {
    OTAssign< TinyMatrix<T,D1,D2> , TinyMatrix<T,D1,D2> ,OpAssign>
    ::apply(*this,rhs, OpAssign());
  }

  // Templated TinyMatrix constructor, maybe to be removed
  template<class T1, unsigned D21, unsigned D22>
  TinyMatrix(const TinyMatrix<T1,D21,D22> &rhs)
  {
    for (unsigned d=0; d<D1*D2; ++d)
      X[d] = (d < D12*D22) ? rhs[d] : T1(0);
  }

  // Constructor from a single T
  TinyMatrix(const T& x00)
  {
    OTAssign<TinyMatrix<T,D1,D2>,T,OpAssign>::apply(*this,x00,OpAssign());
  }

  TinyMatrix(T* x)
  {
    for(int i=0; i<D1*D2; i++)
      X[i] = x[i];
  }

  // Constructors for fixed dimension
  // 1x2 or 2x1
  TinyMatrix(const T& x00, const T& x01)
  {
    X[0] = x00;
    X[1] = x01;
  }

  // 1x3 or 3x1
  TinyMatrix(const T& x00, const T& x01, const T& x02)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
  }

  // 1x4 or 4x1 or 2x2
  TinyMatrix(const T& x00, const T& x01, const T& x02, const T& x03)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
    X[3] = x03;
  }

  // 2x4 or 4x2
  TinyMatrix(const T& x00, const T& x01, const T& x02, const T& x03,
             const T& x04, const T& x05, const T& x06, const T& x07)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
    X[3] = x03;
    X[4] = x04;
    X[5] = x05;
    X[6] = x06;
    X[7] = x07;
  }


  // 3x3
  TinyMatrix(const T& x00, const T& x01, const T& x02, const T& x03,
             const T& x04, const T& x05, const T& x06, const T& x07,
             const T& x08)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
    X[3] = x03;
    X[4] = x04;
    X[5] = x05;
    X[6] = x06;
    X[7] = x07;
    X[8] = x08;
  }

  // 4x4
  TinyMatrix(const T& x00, const T& x01, const T& x02, const T& x03,
             const T& x04, const T& x05, const T& x06, const T& x07,
             const T& x08, const T& x09, const T& x10, const T& x11,
             const T& x12, const T& x13, const T& x14, const T& x15)
  {
    X[ 0] = x00;
    X[ 1] = x01;
    X[ 2] = x02;
    X[ 3] = x03;
    X[ 4] = x04;
    X[ 5] = x05;
    X[ 6] = x06;
    X[ 7] = x07;
    X[ 8] = x08;
    X[ 9] = x09;
    X[10] = x10;
    X[11] = x11;
    X[12] = x12;
    X[13] = x13;
    X[14] = x14;
    X[15] = x15;
  }


  // D1,D2estructor
  ~TinyMatrix() { }

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

  inline TinyMatrix<T,D1,D2>& operator=(const TinyMatrix<T,D1,D2> &rhs)
  {
    OTAssign<TinyMatrix<T,D1,D2>,TinyMatrix<T,D1,D2>,OpAssign>::apply(*this,rhs,OpAssign());
    return *this;
  }

  template<class T1>
  inline TinyMatrix<T,D1,D2>& operator=(const TinyMatrix<T1,D1,D2> &rhs)
  {
    OTAssign<TinyMatrix<T,D1,D2>,TinyMatrix<T1,D1,D2>,OpAssign>::apply(*this,rhs,OpAssign());
    return *this;
  }

  inline TinyMatrix<T,D1,D2>& operator=(const T& rhs)
  {
    OTAssign<TinyMatrix<T,D1,D2>,T,OpAssign>::apply(*this, rhs, OpAssign());
    return *this;
  }

  // Get and Set Operations
  inline Type_t& operator[](unsigned int i)
  {
    return X[i];
  }
  inline Type_t operator[](unsigned int i) const
  {
    return X[i];
  }
  inline Type_t& operator()(unsigned int i)
  {
    return X[i];
  }
  inline Type_t operator()( unsigned int i) const
  {
    return X[i];
  }

  inline Type_t operator()( unsigned int i,  unsigned int j ) const
  {
    return X[i*D2+j];
  }

  inline Type_t& operator()( unsigned int i, unsigned int j )
  {
    return X[i*D2+j];
  }

  // Comparison operators.
  //bool operator==(const TinyMatrix<T,D1,D2>& that) const {
  //  return MetaCompareArrays<T,T,D1,D2>::apply(X,that.X);
  //}
  //bool operator!=(const TinyMatrix<T,D1,D2>& that) const {
  //  return !(*this == that);
  //}

  //----------------------------------------------------------------------
  // parallel communication

  //Message& putMessage(Message& m) const {
  //  m.setCopy(true);
  //  ::putMessage(m, X, X + D1,D2);
  //    return m;
  //}

  //Message& getMessage(Message& m) {
  //  ::getMessage(m, X, X + D1,D2);
  //  return m;
  //}

private:

  // Just store D1,D2 elements of type T.
  T X[D1*D2];
};

// Adding binary operators using macro defined in OhmmsTinyMeta.h
OHMMS_TINYMAT_ACCUM_OPERATORS(TinyMatrix,operator+=,OpAddAssign)
OHMMS_TINYMAT_ACCUM_OPERATORS(TinyMatrix,operator-=,OpSubtractAssign)

OHMMS_TINYMAT_BINARY_OPERATORS(TinyMatrix,operator+,OpAdd)
OHMMS_TINYMAT_BINARY_OPERATORS(TinyMatrix,operator-,OpSubtract)

//----------------------------------------------------------------------
// matrix-matrix multiplication
//----------------------------------------------------------------------
template < class T1, class T2, unsigned D1, unsigned D2, unsigned D3>
inline TinyMatrix<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D1, D3>
dot(const TinyMatrix<T1,D1,D2> &lhs, const TinyMatrix<T2,D2,D3> &rhs)
{
  return OTDot< TinyMatrix<T1,D1,D2> , TinyMatrix<T2,D2,D3> > :: apply(lhs,rhs);
}

template < class T1, class T2, unsigned D1, unsigned D2, unsigned D3>
inline TinyMatrix<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D1, D3>
operator*(const TinyMatrix<T1,D1,D2> &lhs, const TinyMatrix<T2,D2,D3> &rhs)
{
  return OTDot< TinyMatrix<T1,D1,D2> , TinyMatrix<T2,D2,D3> > :: apply(lhs,rhs);
}

//----------------------------------------------------------------------
// matrix-vector multiplication
//----------------------------------------------------------------------
template < class T1, class T2, unsigned D1, unsigned D2>
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D1>
operator*(const TinyMatrix<T1,D1,D2> &lhs, const TinyVector<T2,D2> &rhs)
{
  return OTDot< TinyMatrix<T1,D1,D2> , TinyVector<T2,D2> > :: apply(lhs,rhs);
}


template < class T1, class T2, unsigned D1, unsigned D2>
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D2>
operator*(const TinyVector<T1,D1> &lhs, const TinyMatrix<T2,D1,D2> &rhs)
{
  return OTDot< TinyVector<T1,D1>, TinyMatrix<T2,D1,D2> > :: apply(lhs,rhs);
}

//----------------------------------------------------------------------
// OpMultiply(Matrix,scalar) = element by element operators
// OpMultiply(scalar,Matrix) = element by element operators
// Using OTBinary defined in OhmmsTinyMeta.h
//----------------------------------------------------------------------
template < class T1, class T2, unsigned D1, unsigned D2>
inline TinyMatrix<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D1, D2>
operator*(const TinyMatrix<T1,D1,D2> &lhs, T2 rhs)
{
  return OTBinary< TinyMatrix<T1,D1,D2> , T2, OpMultiply >
         ::apply(lhs,rhs, OpMultiply());
}

template < class T1, class T2, unsigned D1, unsigned D2>
inline TinyMatrix<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D1, D2>
operator*(T1 lhs, const TinyMatrix<T2,D1,D2> &rhs)
{
  return OTBinary< T1, TinyMatrix<T2,D1,D2>, OpMultiply >
         :: apply(lhs,rhs, OpMultiply());
}

//----------------------------------------------------------------------
// I/O
template<class T, unsigned D1, unsigned D2>
ostream& operator<<(std::ostream& out, const TinyMatrix<T,D1,D2>& rhs)
{
  int ii=0;
  for(int i=0; i<D1; i++)
  {
    for(int j=0; j<D2; j++)
      out << rhs[ii++] << " ";
    out << std::endl;
  }
  return out;
}
#endif // OHMMS_TINYMATRIX_H

