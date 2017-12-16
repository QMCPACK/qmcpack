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


#ifndef	SYM_TENZOR_H
#define	SYM_TENZOR_H
/***************************************************************************
 *
 * The POOMA Framework
 *
 * This program was prepared by the Regents of the University of
 * California at Los Alamos National Laboratory (the University) under
 * Contract No.  W-7405-ENG-36 with the U.S. Department of Energy (DOE).
 * The University has certain rights in the program pursuant to the
 * contract and the program should not be copied or distributed outside
 * your organization.  All rights in the program are reserved by the DOE
 * and the University.  Neither the U.S.  Government nor the University
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMA for more details
 *
 ***************************************************************************/



// include files
//#include "Message/Message.h"
#include "PETE/PETE.h"
#include "OhmmsPETE/OhmmsTinyMeta.h"
#include "OhmmsPETE/Tensor.h"

#include <iostream>

namespace qmcplusplus
{

//////////////////////////////////////////////////////////////////////
//
// Definition of class SymTensor.
//
//////////////////////////////////////////////////////////////////////

//
//		| xOO x10 x20 |
//		| x10 x11 x21 |
//		| x20 x21 x22 |
//

template<class T, unsigned D>
class SymTensor
{
public:

  typedef T Type_t;
  enum { ElemDim = 2 };
  enum { Size = D*(D+1)/2 };

  // Default Constructor
  SymTensor()
  {
    OTAssign<SymTensor<T,D>,T,OpAssign>::apply(*this,T(0),OpAssign());
  }

  // A noninitializing ctor.
  class DontInitialize {};
  SymTensor(DontInitialize) {}

  // construct a SymTensor from a single T
  SymTensor(const T& x00)
  {
    OTAssign< SymTensor<T,D>,T,OpAssign>::apply(*this,x00, OpAssign());
  }

  // construct a 2D SymTensor
  SymTensor(const T& x00, const T& x10, const T& x11)
  {
    X[0] = x00;
    X[1] = x10;
    X[2] = x11;
  }
  // construct a 3D SymTensor
  SymTensor(const T& x00, const T& x10, const T& x11, const T& x20,
            const T& x21, const T& x22)
  {
    X[0]= x00;
    X[1]= x10;
    X[2]= x11;
    X[3]= x20;
    X[4]= x21;
    X[5]= x22;
  }

  // Copy Constructor
  template<typename T1>
    SymTensor(const SymTensor<T1,D>& rhs)
    {
      OTAssign< SymTensor<T,D> , SymTensor<T1,D> ,OpAssign > ::
        apply(*this,rhs, OpAssign());
    }

  // Construct from a Tensor.
  // Extract the symmetric part.
  SymTensor(const Tensor<T,D>& t)
  {
    for (int i=0; i<D; ++i)
    {
      (*this)(i,i) = t(i,i);
      for (int j=i+1; j<D; ++j)
        (*this)(i,j) = (t(i,j)+t(j,i))*0.5;
    }
  }

  // Dtor doesn't need to do anything.
  ~SymTensor() { };

  // assignment operators
  const SymTensor<T,D>& operator= (const SymTensor<T,D> &rhs)
  {
    OTAssign< SymTensor<T,D> , SymTensor<T,D> ,OpAssign >
    ::apply(*this,rhs, OpAssign());
    return *this;
  }
  template<class T1>
  const SymTensor<T,D>& operator= (const SymTensor<T1,D> &rhs)
  {
    OTAssign< SymTensor<T,D> , SymTensor<T1,D> ,OpAssign > ::
    apply(*this,rhs, OpAssign());
    return *this;
  }
  const SymTensor<T,D>& operator= (const T& rhs)
  {
    OTAssign< SymTensor<T,D> , T ,OpAssign > :: apply(*this,rhs, OpAssign());
    return *this;
  }
  const SymTensor<T,D>& operator= (const Tensor<T,D> &rhs)
  {
    for (int i=0; i<D; ++i)
    {
      (*this)(i,i) = rhs(i,i);
      for (int j=i+1; j<D; ++j)
        (*this)(i,j) = (rhs(i,j)+rhs(j,i))*0.5;
    }
    return *this;
  }

  // accumulation operators
  template<class T1>
  SymTensor<T,D>& operator+=(const SymTensor<T1,D> &rhs)
  {
    OTAssign< SymTensor<T,D> , SymTensor<T1,D> , OpAddAssign > ::
    apply(*this,rhs, OpAddAssign());
    return *this;
  }
  SymTensor<T,D>& operator+=(const T& rhs)
  {
    OTAssign< SymTensor<T,D>,T,OpAddAssign >::apply(*this,rhs, OpAddAssign());
    return *this;
  }

  template<class T1>
  SymTensor<T,D>& operator-=(const SymTensor<T1,D> &rhs)
  {
    OTAssign< SymTensor<T,D> , SymTensor<T1,D> , OpSubtractAssign > ::
    apply(*this,rhs, OpSubtractAssign());
    return *this;
  }
  SymTensor<T,D>& operator-=(const T& rhs)
  {
    OTAssign< SymTensor<T,D> , T , OpSubtractAssign > ::
    apply(*this,rhs, OpSubtractAssign());
    return *this;
  }

  template<class T1>
  SymTensor<T,D>& operator*=(const SymTensor<T1,D> &rhs)
  {
    OTAssign< SymTensor<T,D> , SymTensor<T1,D> , OpMultiplyAssign > ::
    apply(*this,rhs, OpMultiplyAssign());
    return *this;
  }
  SymTensor<T,D>& operator*=(const T& rhs)
  {
    OTAssign< SymTensor<T,D> , T , OpMultiplyAssign > ::
    apply(*this,rhs, OpMultiplyAssign());
    return *this;
  }

  template<class T1>
  SymTensor<T,D>& operator/=(const SymTensor<T1,D> &rhs)
  {
    OTAssign< SymTensor<T,D> , SymTensor<T1,D> , OpDivideAssign > ::
    apply(*this,rhs, OpDivideAssign());
    return *this;
  }
  SymTensor<T,D>& operator/=(const T& rhs)
  {
    OTAssign< SymTensor<T,D> , T , OpDivideAssign > ::
    apply(*this,rhs, OpDivideAssign());
    return *this;
  }

  // Methods
  void diagonal(const T& rhs)
  {
    for ( int i = 0 ; i < D ; i++ )
    {
      X[((i+1)*i/2) + i] = rhs;
    }
  }

  int len(void)  const
  {
    return Size;
  }
  int size(void) const
  {
    return sizeof(*this);
  }
  int get_Size(void) const
  {
    return Size;
  }

  // Operators

  Type_t operator()(unsigned int i, unsigned int j) const
  {
    int lo = i < j ? i : j;
    int hi = i > j ? i : j;
    return X[((hi+1)*hi/2) + lo];
  }

  Type_t& operator()(unsigned int i, unsigned int j)
  {
    int lo = i < j ? i : j;
    int hi = i > j ? i : j;
    return X[((hi+1)*hi/2) + lo];
  }

//    Type_t& operator()(std::pair<int,int> a) {
//      int i = a.first;
//      int j = a.second;
//      int lo = i < j ? i : j;
//      int hi = i > j ? i : j;
//      return X[((hi+1)*hi/2) + lo];
//    }

//    Type_t operator()( std::pair<int,int> a) const {
//      int i = a.first;
//      int j = a.second;
//      int lo = i < j ? i : j;
//      int hi = i > j ? i : j;
//      return X[((hi+1)*hi/2) + lo];
//    }

  Type_t HL(unsigned int hi, unsigned int lo) const
  {
    PAssert( hi >= lo );
    PAssert( hi<D );
    PAssert( lo<D );
    return X[hi*(hi+1)/2 + lo];
  }
  Type_t& HL(unsigned int hi, unsigned int lo)
  {
    PAssert( hi >= lo );
    PAssert( hi<D );
    PAssert( lo<D );
    return X[hi*(hi+1)/2 + lo];
  }

  Type_t& operator[](unsigned int i)
  {
    PAssert (i < Size);
    return X[i];
  }

  Type_t operator[](unsigned int i) const
  {
    PAssert (i < Size);
    return X[i];
  }

  //TJW: add these 12/16/97 to help with NegReflectAndZeroFace BC:
  // These are the same as operator[] but with () instead:

  Type_t& operator()(unsigned int i)
  {
    PAssert (i < Size);
    return X[i];
  }

  Type_t operator()(unsigned int i) const
  {
    PAssert (i < Size);
    return X[i];
  }
  //TJW.

//    //----------------------------------------------------------------------
//    // Comparison operators.
//    bool operator==(const SymTensor<T,D>& that) const {
//      return MetaCompareArrays<T,T,D*(D+1)/2>::apply(X,that.X);
//    }
//    bool operator!=(const SymTensor<T,D>& that) const {
//      return !(*this == that);
//    }

//    //----------------------------------------------------------------------
//    // parallel communication
//    Message& putMessage(Message& m) const {
//      m.setCopy(true);
//      ::putMessage(m, X, X + ((D*(D + 1)/2)));
//      return m;
//    }

//    Message& getMessage(Message& m) {
//      ::getMessage(m, X, X + ((D*(D + 1)/2)));
//      return m;
//    }

private:

  // The elements themselves.
  T X[Size];

};


//////////////////////////////////////////////////////////////////////
//
// Free functions
//
//////////////////////////////////////////////////////////////////////

template <class T, unsigned D>
T trace(const SymTensor<T,D> &rhs)
{
  T result = 0.0;
  for (int i = 0 ; i < D ; i++ )
    result += rhs(i,i);
  return result;
}

template <class T, unsigned D>
SymTensor<T,D> transpose(const SymTensor<T,D> &rhs)
{
  return rhs;
}

//////////////////////////////////////////////////////////////////////
//
// Unary Operators
//
//////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------
// unary operator-
//template<class T, unsigned D>
//inline SymTensor<T,D> operator-(const SymTensor<T,D> &op)
//{
//  return MetaUnary< SymTensor<T,D> , OpUnaryMinus > :: apply(op,OpUnaryMinus());
//}

//----------------------------------------------------------------------
// unary operator+
//template<class T, unsigned D>
//inline const SymTensor<T,D> &operator+(const SymTensor<T,D> &op)
//{
//  return op;
//}

//////////////////////////////////////////////////////////////////////
//
// Binary Operators
//
//////////////////////////////////////////////////////////////////////

//
// Elementwise operators.
//

OHMMS_META_BINARY_OPERATORS(SymTensor,operator+,OpAdd)
OHMMS_META_BINARY_OPERATORS(SymTensor,operator-,OpSubtract)
OHMMS_META_BINARY_OPERATORS(SymTensor,operator*,OpMultiply)
OHMMS_META_BINARY_OPERATORS(SymTensor,operator/,OpDivide)

/*
TSV_ELEMENTWISE_OPERATOR2(SymTensor,Tensor,operator+,OpAdd)
TSV_ELEMENTWISE_OPERATOR2(Tensor,SymTensor,operator+,OpAdd)
TSV_ELEMENTWISE_OPERATOR2(SymTensor,Tensor,operator-,OpSubtract)
TSV_ELEMENTWISE_OPERATOR2(Tensor,SymTensor,operator-,OpSubtract)
TSV_ELEMENTWISE_OPERATOR2(SymTensor,Tensor,operator*,OpMultiply)
TSV_ELEMENTWISE_OPERATOR2(Tensor,SymTensor,operator*,OpMultiply)
TSV_ELEMENTWISE_OPERATOR2(SymTensor,Tensor,operator/,OpDivide)
TSV_ELEMENTWISE_OPERATOR2(Tensor,SymTensor,operator/,OpDivide)
*/

//----------------------------------------------------------------------
// dot products.
//----------------------------------------------------------------------

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const SymTensor<T1,D> &lhs, const SymTensor<T2,D> &rhs)
{
  return OTDot< SymTensor<T1,D> , SymTensor<T2,D> > :: apply(lhs,rhs);
}

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const SymTensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
{
  return OTDot< SymTensor<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
}

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const Tensor<T1,D> &lhs, const SymTensor<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , SymTensor<T2,D> > :: apply(lhs,rhs);
}

template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const TinyVector<T1,D> &lhs, const SymTensor<T2,D> &rhs)
{
  return OTDot< TinyVector<T1,D> , SymTensor<T2,D> > :: apply(lhs,rhs);
}

template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const SymTensor<T1,D> &lhs, const TinyVector<T2,D> &rhs)
{
  return OTDot< SymTensor<T1,D> , TinyVector<T2,D> > :: apply(lhs,rhs);
}

//----------------------------------------------------------------------
// double dot products.
//----------------------------------------------------------------------

//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const SymTensor<T1,D> &lhs, const SymTensor<T2,D> &rhs)
//{
//  return MetaDotDot< SymTensor<T1,D> , SymTensor<T2,D> > :: apply(lhs,rhs);
//}
//
//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const SymTensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
//{
//  return MetaDotDot< SymTensor<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
//}
//
//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const Tensor<T1,D> &lhs, const SymTensor<T2,D> &rhs)
//{
//  return MetaDotDot< Tensor<T1,D> , SymTensor<T2,D> > :: apply(lhs,rhs);
//}
//----------------------------------------------------------------------
// I/O
template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const SymTensor<T,D>& rhs)
{
  if (D >= 1)
  {
    for (int i=0; i<D; i++)
    {
      for (int j=0; j<D-1; j++)
      {
        out << rhs(i,j) << " ";
      }
      out << rhs(i,D-1) << " ";
      if (i < D - 1)
        out << std::endl;
    }
  }
  else
  {
    out << " " << rhs(0,0) << " ";
  }
  return out;
}

}
#endif // SYM_TENZOR_H


