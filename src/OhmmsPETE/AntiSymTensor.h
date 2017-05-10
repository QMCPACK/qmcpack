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


#ifndef	ANTI_SYM_TENZOR_H
#define	ANTI_SYM_TENZOR_H

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
//#include "Utility/PAssert.h"
//#include "Message/Message.h"
#include "PETE/PETE.h"
#include "OhmmsPETE/OhmmsTinyMeta.h"
#include "OhmmsPETE/Tensor.h"

#include <iostream>

namespace qmcplusplus
{

//////////////////////////////////////////////////////////////////////
//
// Definition of class AntiSymTensor.
//
//////////////////////////////////////////////////////////////////////

//
//		|  O -x10  -x20 |
//		| x10   0  -x21 |
//		| x20  x21   0  |
//

template<class T, unsigned D>
class AntiSymTensor
{
public:

  typedef T Type_t;
  enum { ElemDim = 2 };
  enum { Size = D*(D-1)/2 };

  // Default Constructor
  AntiSymTensor()
  {
    OTAssign< AntiSymTensor<T,D>,T,OpAssign>::apply(*this,T(0),OpAssign());
  }

  // A noninitializing ctor.
  class DontInitialize {};
  AntiSymTensor(DontInitialize) {}

  // Construct an AntiSymTensor from a single T.
  // This doubles as the 2D AntiSymTensor initializer.
  AntiSymTensor(const T& x00)
  {
    OTAssign< AntiSymTensor<T,D>,T,OpAssign>::apply(*this,x00,OpAssign());
  }
  // construct a 3D AntiSymTensor
  AntiSymTensor(const T& x10, const T& x20, const T& x21)
  {
    X[0]= x10;
    X[1]= x20;
    X[2]= x21;
  }

  // Copy Constructor
  AntiSymTensor(const AntiSymTensor<T,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T,D> ,OpAssign > ::
    apply(*this,rhs,OpAssign());
  }

  // Construct from a Tensor.
  // Extract the antisymmetric part.
  AntiSymTensor( const Tensor<T,D>& t )
  {
    for (int i=1; i<D; ++i)
    {
      for (int j=0; j<i; ++j)
        (*this)[((i-1)*i/2)+j] = (t(i,j)-t(j,i))*0.5;
    }
  }

  ~AntiSymTensor() { };

  // assignment operators
  const AntiSymTensor<T,D>& operator= (const AntiSymTensor<T,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T,D> ,OpAssign > ::
    apply(*this,rhs,OpAssign());
    return *this;
  }
  template<class T1>
  const AntiSymTensor<T,D>& operator= (const AntiSymTensor<T1,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T1,D> ,OpAssign > ::
    apply(*this,rhs,OpAssign());
    return *this;
  }
  const AntiSymTensor<T,D>& operator= (const T& rhs)
  {
    OTAssign< AntiSymTensor<T,D> , T ,OpAssign > ::
    apply(*this,rhs,OpAssign());
    return *this;
  }

  // accumulation operators
  template<class T1>
  AntiSymTensor<T,D>& operator+=(const AntiSymTensor<T1,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T1,D> , OpAddAssign >
    :: apply(*this,rhs,OpAddAssign());
    return *this;
  }

  template<class T1>
  AntiSymTensor<T,D>& operator-=(const AntiSymTensor<T1,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T1,D> ,
              OpSubtractAssign > :: apply(*this,rhs,OpSubtractAssign());
    return *this;
  }

  template<class T1>
  AntiSymTensor<T,D>& operator*=(const AntiSymTensor<T1,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T1,D> ,
              OpMultiplyAssign > :: apply(*this,rhs,OpMultiplyAssign());
    return *this;
  }
  AntiSymTensor<T,D>& operator*=(const T& rhs)
  {
    OTAssign< AntiSymTensor<T,D> , T , OpMultiplyAssign > ::
    apply(*this,rhs,OpMultiplyAssign());
    return *this;
  }

  template<class T1>
  AntiSymTensor<T,D>& operator/=(const AntiSymTensor<T1,D> &rhs)
  {
    OTAssign< AntiSymTensor<T,D> , AntiSymTensor<T1,D> ,
              OpDivideAssign > :: apply(*this,rhs,OpDivideAssign());
    return *this;
  }
  AntiSymTensor<T,D>& operator/=(const T& rhs)
  {
    OTAssign< AntiSymTensor<T,D> , T , OpDivideAssign > ::
    apply(*this,rhs,OpDivideAssign());
    return *this;
  }

  // Methods

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

  class AssignProxy
  {
  public:
    AssignProxy(Type_t &elem, int where)
      : elem_m(elem), where_m(where) { }
    AssignProxy(const AssignProxy &model)
      : elem_m(model.elem_m), where_m(where_m) { }
    const AssignProxy &operator=(const AssignProxy &a)
    {
      PAssert(where_m != 0 || a.elem_m == -a.elem_m);
      elem_m = where_m < 0 ? -a.elem_m : a.elem_m;
      return *this;
    }
    const AssignProxy &operator=(const Type_t &e)
    {
      PAssert(where_m != 0 || e == -e);
      elem_m = where_m < 0 ? -e : e;
      return *this;
    }

    operator Type_t() const
    {
      return (where_m < 0 ? -elem_m : elem_m);
    }

  private:

    //mutable Type_t &elem_m;
    //mutable int where_m;
    Type_t &elem_m;
    int where_m;
  };

  // Operators

  Type_t operator()(unsigned int i, unsigned int j) const
  {
    if (i == j)
      return T(0.0);
    else
      if (i < j)
        return -X[((j-1)*j/2) + i];
      else
        return X[((i-1)*i/2) + j];
  }

//    Type_t operator()( std::pair<int,int> a) const {
//      int i = a.first;
//      int j = a.second;
//      return (*this)(i, j);
//    }

  AssignProxy operator()(unsigned int i, unsigned int j)
  {
    if (i == j)
      return AssignProxy(AntiSymTensor<T,D>::Zero, 0);
    else
    {
      int lo = i < j ? i : j;
      int hi = i > j ? i : j;
      return AssignProxy(X[((hi-1)*hi/2) + lo], i - j);
    }
  }

//    AssignProxy operator()(std::pair<int,int> a) {
//      int i = a.first;
//      int j = a.second;
//      return (*this)(i, j);
//    }

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

//    //----------------------------------------------------------------------
//    // Comparison operators.
//    bool operator==(const AntiSymTensor<T,D>& that) const {
//      return TSV_MetaCompareArrays<T,T,D*(D-1)/2>::apply(X,that.X);
//    }
//    bool operator!=(const AntiSymTensor<T,D>& that) const {
//      return !(*this == that);
//    }

//    //----------------------------------------------------------------------
//    // parallel communication
//    Message& putMessage(Message& m) const {
//      m.setCopy(true);
//      ::putMessage(m, X, X + Size);
//      return m;
//    }

//    Message& getMessage(Message& m) {
//      ::getMessage(m, X, X + Size);
//      return m;
//    }

private:

  friend class AssignProxy;

  // The elements themselves.
  T X[Size];

  // A place to store a zero element.
  // We need to return a reference to this
  // for the diagonal element.
  static T Zero;
};


///////////////////////////////////////////////////////////////////////////
// Specialization for 1D  -- this is basically just the zero tensor
///////////////////////////////////////////////////////////////////////////

template <class T>
class AntiSymTensor<T,1>
{
public:

  typedef T Type_t;
  enum { ElemDim = 2 };

  // Default Constructor
  AntiSymTensor() {}

  // Copy Constructor
  AntiSymTensor(const AntiSymTensor<T,1>&) {}

  ~AntiSymTensor() {}

  // assignment operators
  const AntiSymTensor<T,1>& operator=(const AntiSymTensor<T,1>&)
  {
    return *this;
  }
  template <class T1>
  const AntiSymTensor<T,1>& operator=(const AntiSymTensor<T1,1>&)
  {
    return *this;
  }
  const AntiSymTensor<T,1>& operator=(const T& rhs)
  {
    return *this;
  }

  // accumulation operators
  template <class T1>
  AntiSymTensor<T,1>& operator+=(const AntiSymTensor<T1,1>&)
  {
    return *this;
  }

  template <class T1>
  AntiSymTensor<T,1>& operator-=(const AntiSymTensor<T1,1>&)
  {
    return *this;
  }

  template <class T1>
  AntiSymTensor<T,1>& operator*=(const AntiSymTensor<T1,1>&)
  {
    return *this;
  }
  AntiSymTensor<T,1>& operator*=(const T&)
  {
    return *this;
  }

  template <class T1>
  AntiSymTensor<T,1>& operator/=(const AntiSymTensor<T1,1>&)
  {
    return *this;
  }
  AntiSymTensor<T,1>& operator/=(const T&)
  {
    return *this;
  }

  // Methods

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

  class AssignProxy
  {
  public:
    AssignProxy(Type_t &elem, int where)
      : elem_m(elem), where_m(where) {}
    AssignProxy(const AssignProxy& model)
      : elem_m(model.elem_m), where_m(where_m) {}
    const AssignProxy& operator=(const AssignProxy& a)
    {
      PAssert(where_m != 0 || a.elem_m == -a.elem_m);
      elem_m = where_m < 0 ? -a.elem_m : a.elem_m;
      return *this;
    }
    const AssignProxy& operator=(const Type_t& e)
    {
      PAssert(where_m != 0 || e == -e);
      elem_m = where_m < 0 ? -e : e;
      return *this;
    }

    operator Type_t() const
    {
      return (where_m < 0 ? -elem_m : elem_m);
    }

  private:

    //mutable Type_t &elem_m;
    //mutable int where_m;
    Type_t &elem_m;
    int where_m;
  };

  // Operators

  Type_t operator()(unsigned int i, unsigned int j) const
  {
    PAssert(i==j);
    return T(0.0);
  }

//    Type_t operator()( std::pair<int,int> a) const {
//      int i = a.first;
//      int j = a.second;
//      return (*this)(i, j);
//    }

  AssignProxy operator()(unsigned int i, unsigned int j)
  {
    PAssert(i==j);
    return AssignProxy(AntiSymTensor<T,1>::Zero, 0);
  }

//    AssignProxy operator()(std::pair<int,int> a) {
//      int i = a.first;
//      int j = a.second;
//      return (*this)(i, j);
//    }

  Type_t operator[](unsigned int i) const
  {
    PAssert (i == 0);
    return AntiSymTensor<T,1>::Zero;
  }

  // These are the same as operator[] but with () instead:

  Type_t operator()(unsigned int i) const
  {
    PAssert (i == 0);
    return AntiSymTensor<T,1>::Zero;
  }

  //----------------------------------------------------------------------
  // Comparison operators.
  bool operator==(const AntiSymTensor<T,1>& that) const
  {
    return true;
  }
  bool operator!=(const AntiSymTensor<T,1>& that) const
  {
    return !(*this == that);
  }

  //----------------------------------------------------------------------
  // parallel communication
//    Message& putMessage(Message& m) const {
//      m.setCopy(true);
//      m.put(AntiSymTensor<T,1>::Zero);
//      return m;
//    }

//    Message& getMessage(Message& m) {
//      T zero;
//      m.get(zero);
//      return m;
//    }

private:

  friend class AssignProxy;

  // The number of elements.
  enum { Size = 0 };

  // A place to store a zero element.
  // We need to return a reference to this
  // for the diagonal element.
  static T Zero;
};


//////////////////////////////////////////////////////////////////////
//
// Free functions
//
//////////////////////////////////////////////////////////////////////

template <class T, unsigned D>
T trace(const AntiSymTensor<T,D>&)
{
  return T(0.0);
}

template <class T, unsigned D>
AntiSymTensor<T,D> transpose(const AntiSymTensor<T,D>& rhs)
{
  return -rhs;
}
//////////////////////////////////////////////////////////////////////
//
// Unary Operators
//
//////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------
// unary operator-
//  template<class T, unsigned D>
//  inline AntiSymTensor<T,D> operator-(const AntiSymTensor<T,D> &op)
//  {
//    return MetaUnary< AntiSymTensor<T,D> , OpUnaryMinus > :: apply(op,OpUnaryMinus());
//  }

//----------------------------------------------------------------------
// unary operator+
template<class T, unsigned D>
inline const AntiSymTensor<T,D> &operator+(const AntiSymTensor<T,D> &op)
{
  return op;
}

//////////////////////////////////////////////////////////////////////
//
// Binary Operators
//
//////////////////////////////////////////////////////////////////////

//
// Elementwise operators.
//

OHMMS_META_BINARY_OPERATORS(AntiSymTensor,operator+,OpAdd)
OHMMS_META_BINARY_OPERATORS(AntiSymTensor,operator-,OpSubtract)
OHMMS_META_BINARY_OPERATORS(AntiSymTensor,operator*,OpMultiply)
OHMMS_META_BINARY_OPERATORS(AntiSymTensor,operator/,OpDivide)

//----------------------------------------------------------------------
// dot products.
//----------------------------------------------------------------------

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const AntiSymTensor<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
{
  return OTDot< AntiSymTensor<T1,D> , AntiSymTensor<T2,D> > ::
         apply(lhs,rhs);
}

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const AntiSymTensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , Tensor<T2,D> > ::
         apply(Tensor<T1,D>(lhs),rhs);
}

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const Tensor<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , Tensor<T2,D> > ::
         apply(lhs,Tensor<T2,D>(rhs));
}

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const AntiSymTensor<T1,D> &lhs, const SymTensor<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , Tensor<T2,D> > ::
         apply(Tensor<T1,D>(lhs),Tensor<T2,D>(rhs));
}

template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const SymTensor<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , Tensor<T2,D> > ::
         apply(Tensor<T1,D>(lhs),Tensor<T2,D>(rhs));
}

template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const TinyVector<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
{
  return OTDot< TinyVector<T1,D> , AntiSymTensor<T2,D> > :: apply(lhs,rhs);
}

template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const AntiSymTensor<T1,D> &lhs, const TinyVector<T2,D> &rhs)
{
  return OTDot< AntiSymTensor<T1,D> , TinyVector<T2,D> > :: apply(lhs,rhs);
}

//----------------------------------------------------------------------
// double dot products.
//----------------------------------------------------------------------

//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const AntiSymTensor<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
//{
//  return MetaDotDot< AntiSymTensor<T1,D> , AntiSymTensor<T2,D> > ::
//    apply(lhs,rhs);
//}

//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const AntiSymTensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
//{
//  return MetaDotDot< Tensor<T1,D> , Tensor<T2,D> > ::
//    apply(Tensor<T1,D>(lhs),rhs);
//}

//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const Tensor<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
//{
//  return MetaDotDot< Tensor<T1,D> , Tensor<T2,D> > ::
//    apply(lhs,Tensor<T2,D>(rhs));
//}

//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const AntiSymTensor<T1,D> &lhs, const SymTensor<T2,D> &rhs)
//{
//  return MetaDotDot< Tensor<T1,D> , Tensor<T2,D> > ::
//    apply(Tensor<T1,D>(lhs),Tensor<T2,D>(rhs));
//}

//template < class T1, class T2, unsigned D >
//inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
//dotdot(const SymTensor<T1,D> &lhs, const AntiSymTensor<T2,D> &rhs)
//{
//  return MetaDotDot< Tensor<T1,D> , Tensor<T2,D> > ::
//    apply(Tensor<T1,D>(lhs),Tensor<T2,D>(rhs));
//}

//----------------------------------------------------------------------
// I/O
template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const AntiSymTensor<T,D>& rhs)
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

//////////////////////////////////////////////////////////////////////

template<class T, unsigned int D>
T AntiSymTensor<T,D>::Zero = 0;

}
#endif // ANTI_SYM_TENZOR_H


