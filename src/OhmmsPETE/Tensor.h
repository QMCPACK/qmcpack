//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//		      Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//		      Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
// 		      Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef	OHMMS_TENSOR_H
#define	OHMMS_TENSOR_H

#include "PETE/PETE.h"
#include "OhmmsPETE/OhmmsTinyMeta.h"
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



namespace qmcplusplus
{
// forward declarations
template <class T, unsigned D> class SymTensor;
template <class T, unsigned D> class AntiSymTensor;


/** Tensor<T,D>  class for D by D tensor
 *
 * @tparam T datatype
 * @tparm D dimension
 */
template<class T, unsigned D>
class Tensor
{
public:

  typedef T Type_t;
  enum { ElemDim = 2 };
  enum { Size = D*D };

  // Default Constructor
  inline Tensor()
  {
    for(size_t d=0; d<Size; ++d)
      X[d] = T(0);
  }

  // A noninitializing ctor.
  class DontInitialize {};
  inline Tensor(DontInitialize) {}

  // Copy Constructor
  inline Tensor(const Tensor& rhs) = default;

  template<class T1>
  inline Tensor(const Tensor<T1,D> &rhs)
  {
    for(size_t d=0; d<Size; ++d)
      X[d] = rhs[d];
  }

  // constructor from a single T
  inline Tensor(const T& x00)
  {
    for(size_t d=0; d<Size; ++d)
      X[d] = x00;
  }

  // constructors for fixed dimension
  Tensor(const T& x00, const T& x10, const T& x01, const T& x11)
  {
    X[0] = x00;
    X[1] = x10;
    X[2] = x01;
    X[3] = x11;
  }
  Tensor(const T& x00, const T& x10, const T& x20, const T& x01, const T& x11,
         const T& x21, const T& x02, const T& x12, const T& x22)
  {
    X[0] = x00;
    X[1] = x10;
    X[2] = x20;
    X[3] = x01;
    X[4] = x11;
    X[5] = x21;
    X[6] = x02;
    X[7] = x12;
    X[8] = x22;
  }

  //constructor from SymTensor
  Tensor(const SymTensor<T,D>&);

  // constructor from AntiSymTensor
  Tensor(const AntiSymTensor<T,D>&);

  // destructor
  ~Tensor() { };

  // assignment operators
  inline Tensor& operator=(const Tensor& rhs) = default;

  template<class T1>
  inline Tensor<T,D>& operator= (const Tensor<T1,D> &rhs)
  {
    OTAssign< Tensor<T,D> , Tensor<T1,D> ,OpAssign> :: apply(*this,rhs,OpAssign());
    return *this;
  }
  inline Tensor<T,D>& operator= (const T& rhs)
  {
    OTAssign< Tensor<T,D> , T ,OpAssign> :: apply(*this,rhs, OpAssign());
    return *this;
  }

  // accumulation operators
  template<class T1>
  inline Tensor<T,D>& operator+=(const Tensor<T1,D> &rhs)
  {
    OTAssign< Tensor<T,D> , Tensor<T1,D> , OpAddAssign > ::
    apply(*this,rhs,OpAddAssign());
    return *this;
  }
  inline Tensor<T,D>& operator+=(const T& rhs)
  {
    OTAssign< Tensor<T,D> , T , OpAddAssign > ::
    apply(*this,rhs, OpAddAssign());
    return *this;
  }

  template<class T1>
  inline Tensor<T,D>& operator-=(const Tensor<T1,D> &rhs)
  {
    OTAssign< Tensor<T,D> , Tensor<T1,D> , OpSubtractAssign > ::
    apply(*this,rhs,OpSubtractAssign());
    return *this;
  }

  inline Tensor<T,D>& operator-=(const T& rhs)
  {
    OTAssign< Tensor<T,D> , T , OpSubtractAssign > ::
    apply(*this,rhs,OpSubtractAssign());
    return *this;
  }

  template<class T1>
  inline Tensor<T,D>& operator*=(const Tensor<T1,D> &rhs)
  {
    OTAssign< Tensor<T,D> , Tensor<T1,D> , OpMultiplyAssign > ::
    apply(*this,rhs,OpMultiplyAssign());
    return *this;
  }

  inline Tensor<T,D>& operator*=(const T& rhs)
  {
    OTAssign< Tensor<T,D> , T , OpMultiplyAssign > ::
    apply(*this,rhs, OpMultiplyAssign());
    return *this;
  }

  template<class T1>
  inline Tensor<T,D>& operator/=(const Tensor<T1,D> &rhs)
  {
    OTAssign< Tensor<T,D> , Tensor<T1,D> , OpDivideAssign > ::
    apply(*this,rhs, OpDivideAssign());
    return *this;
  }

  inline Tensor<T,D>& operator/=(const T& rhs)
  {
    OTAssign< Tensor<T,D> , T , OpDivideAssign > ::
    apply(*this,rhs, OpDivideAssign());
    return *this;
  }

  // Methods

  inline void diagonal(const T& rhs)
  {
    for (int i = 0 ; i < D ; i++ )
      (*this)(i,i) = rhs;
  }

  inline void add2diagonal(T rhs)
  {
    for (int i = 0 ; i < D ; i++ )
      (*this)(i,i) += rhs;
  }

  ///return the size
  inline int len()  const
  {
    return Size;
  }
  ///return the size
  inline int size() const
  {
    return Size;
  }

  /** return the i-th value or assign
   * @param i index [0,D*D)
   */
  inline Type_t &operator[]( unsigned int i )
  {
    return X[i];
  }

  /** return the i-th value
   * @param i index [0,D*D)
   */
  inline Type_t operator[]( unsigned int i ) const
  {
    return X[i];
  }

  //TJW: add these 12/16/97 to help with NegReflectAndZeroFace BC:
  // These are the same as operator[] but with () instead:
  inline Type_t& operator()(unsigned int i)
  {
    return X[i];
  }

  inline Type_t operator()(unsigned int i) const
  {
    return X[i];
  }
  //TJW.

  /** return the (i,j) component
   * @param i index [0,D)
   * @param j index [0,D)
   */
  inline Type_t operator()( unsigned int i,  unsigned int j ) const
  {
    return X[i*D+j];
  }

  /** return/assign the (i,j) component
   * @param i index [0,D)
   * @param j index [0,D)
   */
  inline Type_t& operator()( unsigned int i, unsigned int j )
  {
    return X[i*D+j];
  }

  inline TinyVector<T,D> getRow(unsigned int i)
  {
    TinyVector<T,D> res;
    for(int j=0; j<D; j++)
      res[j]=X[i*D+j];
    return res;
  }

  inline TinyVector<T,D> getColumn(unsigned int i)
  {
    TinyVector<T,D> res;
    for(int j=0; j<D; j++)
      res[j]=X[j*D+i];
    return res;
  }

  inline Type_t* data()
  {
    return  X;
  }
  inline const Type_t* data() const
  {
    return  X;
  }
  inline Type_t* begin()
  {
    return  X;
  }
  inline const Type_t* begin() const
  {
    return  X;
  }
  inline Type_t* end()
  {
    return  X+Size;
  }
  inline const Type_t* end() const
  {
    return  X+Size;
  }

//  Removed operator using std::pair
//    Type_t operator()(const std::pair<int,int> i) const {
//      PAssert ( (i.first>=0) && (i.second>=0) && (i.first<D) && (i.second<D) );
//      return (*this)(i.first,i.second);
//    }

//    Type_t& operator()(const std::pair<int,int> i) {
//      PAssert ( (i.first>=0) && (i.second>=0) && (i.first<D) && (i.second<D) );
//      return (*this)(i.first,i.second);
//    }


//----------------------------------------------------------------------
// Comparison operators.
//    bool operator==(const Tensor<T,D>& that) const {
//      return MetaCompareArrays<T,T,D*D>::apply(X,that.X);
//    }
//    bool operator!=(const Tensor<T,D>& that) const {
//      return !(*this == that);
//    }

  //----------------------------------------------------------------------
  // parallel communication
//    Message& putMessage(Message& m) const {
//      m.setCopy(true);
//      const T *p = X;
//      ::putMessage(m, p, p + D*D);
//      return m;
//    }

//    Message& getMessage(Message& m) {
//      T *p = X;
//      ::getMessage(m, p, p + D*D);
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

/** trace \f$ result = \sum_k rhs(k,k)\f$
 * @param rhs a tensor
 */
template <class T, unsigned D>
inline T trace(const Tensor<T,D>& rhs)
{
  T result = 0.0;
  for (int i = 0 ; i < D ; i++ )
    result += rhs(i,i);
  return result;
}

/** transpose a tensor
 */
template <class T, unsigned D>
inline Tensor<T,D> transpose(const Tensor<T,D>& rhs)
{
  Tensor<T,D> result; // = Tensor<T,D>::DontInitialize();
  for (int j = 0 ; j < D ; j++ )
    for (int i = 0 ; i < D ; i++ )
      result(i,j) = rhs(j,i);
  return result;
}

/** Tr(a*b), \f$ \sum_i\sum_j a(i,j)*b(j,i) \f$
 */
template <class T1, class T2, unsigned D>
inline T1 trace(const Tensor<T1,D>& a, const Tensor<T2,D>& b)
{
  T1 result = 0.0;
  for (int i = 0 ; i < D ; i++ )
    for(int j=0; j<D; j++)
      result += a(i,j)*b(j,i);
  return result;
}

/** Tr(a^t *b), \f$ \sum_i\sum_j a(i,j)*b(i,j) \f$
 */
template <class T, unsigned D>
inline T traceAtB(const Tensor<T,D>& a, const Tensor<T,D>& b)
{
  T result = 0.0;
  for (int i = 0 ; i < D*D ; i++ )
    result += a(i)*b(i);
  return result;
}

/** Tr(a^t *b), \f$ \sum_i\sum_j a(i,j)*b(i,j) \f$
 */
template <class T1, class T2, unsigned D>
inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t traceAtB(const Tensor<T1,D>& a, const Tensor<T2,D>& b)
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T;
  T result = 0.0;
  for (int i = 0 ; i < D*D ; i++ )
    result += a(i)*b(i);
  return result;
}

//////////////////////////////////////////////////////////////////////
//
// Unary Operators
//
//////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------
// unary operator-
//template<class T, unsigned D>
//inline Tensor<T,D> operator-(const Tensor<T,D> &op)
//{
//  return MetaUnary< Tensor<T,D> , OpUnaryMinus > :: apply(op,OpUnaryMinus());
//}
//----------------------------------------------------------------------
// unary operator+
//template<class T, unsigned D>
//inline const Tensor<T,D> &operator+(const Tensor<T,D> &op)
//{
//  return op;
//}

/// Binary Operators
OHMMS_META_BINARY_OPERATORS(Tensor,operator+,OpAdd)
OHMMS_META_BINARY_OPERATORS(Tensor,operator-,OpSubtract)
OHMMS_META_BINARY_OPERATORS(Tensor,operator*,OpMultiply)
OHMMS_META_BINARY_OPERATORS(Tensor,operator/,OpDivide)
//TSV_ELEMENTWISE_OPERATOR(Tensor,Min,FnMin)
//TSV_ELEMENTWISE_OPERATOR(Tensor,Max,FnMax)

/** Tensor-Tensor dot product \f$result(i,j)=\sum_k lhs(i,k)*rhs(k,j)\f$
 * @param lhs  a tensor
 * @param rhs  a tensor
 */
template < class T1, class T2, unsigned D >
inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D>
dot(const Tensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
}

/** Vector-Tensor dot product \f$result(i)=\sum_k lhs(k)*rhs(k,i)\f$
 * @param lhs  a vector
 * @param rhs  a tensor
 */
template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t, D>
dot(const TinyVector<T1,D> &lhs, const Tensor<T2,D> &rhs)
{
  return OTDot< TinyVector<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
}

/** Tensor-Vector dot product \f$result(i)=\sum_k lhs(i,k)*rhs(k)\f$
 * @param lhs  a tensor
 * @param rhs  a vector
 */
template < class T1, class T2, unsigned D >
inline TinyVector<typename BinaryReturn<T1,T2,OpMultiply>::Type_t,D>
dot(const Tensor<T1,D> &lhs, const TinyVector<T2,D> &rhs)
{
  return OTDot< Tensor<T1,D> , TinyVector<T2,D> > :: apply(lhs,rhs);
}

//----------------------------------------------------------------------
// double dot product.
//----------------------------------------------------------------------

// template < class T1, class T2, unsigned D >
// inline typename BinaryReturn<T1,T2,OpMultiply>::Type_t
// dotdot(const Tensor<T1,D> &lhs, const Tensor<T2,D> &rhs)
// {
//   return MetaDotDot< Tensor<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
//}


//----------------------------------------------------------------------
// Outer product.
//----------------------------------------------------------------------

///** Vector-vector outter product \f$ result(i,j)=v1(i)*v2(j)\f$
// * @param v1 a vector
// * @param v2 a vector
// */
//  template<class T1, class T2, unsigned int D>
//  inline Tensor<typename BinaryReturn<T1,T2,OpMultiply>::type,D >
//  outerProduct(const TinyVector<T1,D>& v1, const TinyVector<T2,D>& v2)
//  {
//    typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t T0;
////    //#if (defined(POOMA_SGI_CC_721_TYPENAME_BUG) || defined(__MWERKS__))
////    //Tensor<T0,D> ret = Tensor<T0,D>::DontInitialize();
////    //#else
//    Tensor<T0,D> ret = typename Tensor<T0,D>::DontInitialize();
////    //#endif // POOMA_SGI_CC_721_TYPENAME_BUG
//    for (unsigned int i=0; i<D; ++i)
//      for (unsigned int j=0; j<D; ++j)
//        ret(i,j) = v1[i]*v2[j];
//    return ret;
//  }
//
//  template<class T1, unsigned int D>
//  inline Tensor<T1,D >
//  outerProduct(const TinyVector<T1,D>& v1, const TinyVector<T1,D>& v2)
//  {
//    Tensor<T1,D> ret = typename Tensor<T1,D>::DontInitialize();
////    //#endif // POOMA_SGI_CC_721_TYPENAME_BUG
//    for (unsigned int i=0; i<D; ++i)
//      for (unsigned int j=0; j<D; ++j)
//        ret(i,j) = v1[i]*v2[j];
//    return ret;
//  }
//


//----------------------------------------------------------------------
// I/O
template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const Tensor<T,D>& rhs)
{
  if (D >= 1)
  {
    for (int i=0; i<D; i++)
    {
      for (int j=0; j<D-1; j++)
      {
        out << rhs(i,j) << "  ";
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

template<class T, unsigned D>
std::istream& operator>>(std::istream& is, Tensor<T,D>& rhs)
{
  for(int i=0; i<D*D; i++)
    is >> rhs[i];
  return is;
}

}
// include header files for SymTensor and AntiSymTensor in order
// to define constructors for Tensor using these types
#include "OhmmsPETE/SymTensor.h"
#include "OhmmsPETE/AntiSymTensor.h"

namespace qmcplusplus
{
template <class T, unsigned D>
Tensor<T,D>::Tensor(const SymTensor<T,D>& rhs)
{
  for (int i=0; i<D; ++i)
    for (int j=0; j<D; ++j)
      (*this)(i,j) = rhs(i,j);
}

template <class T, unsigned D>
Tensor<T,D>::Tensor(const AntiSymTensor<T,D>& rhs)
{
  for (int i=0; i<D; ++i)
    for (int j=0; j<D; ++j)
      (*this)(i,j) = rhs(i,j);
}

}

#endif // OHMMS_TENSOR_H

