//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_TINYVECTOR_H
#define OHMMS_TINYVECTOR_H

/***************************************************************************
 *
 * \class TinyVector
 * \brief Pooma/AppyTypes/Vecktor is modified to work with PETE.
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
#include <iomanip>
#include "PETE/PETE.h"
#include "OhmmsPETE/OhmmsTinyMeta.h"

namespace qmcplusplus
{
/** Fixed-size array. candidate for array<T,D>
 */
template<class T, unsigned D>
struct TinyVector
{
  using Type_t = T;
  enum
  {
    Size = D
  };
  T X[Size];

  // Default Constructor initializes to zero.
  inline TinyVector()
  {
    for (size_t d = 0; d < D; ++d)
      X[d] = T(0);
  }

  // Copy Constructor
  inline TinyVector(const TinyVector& rhs) = default;

  // Move Constructor
  inline TinyVector(TinyVector&& rhs) = default;

  // Templated TinyVector constructor.
  template<class T1>
  inline TinyVector(const TinyVector<T1, D>& rhs)
  {
    for (size_t d = 0; d < D; ++d)
      X[d] = rhs[d];
  }

  // Constructor from a single T
  inline TinyVector(const T& x00)
  {
    for (size_t d = 0; d < D; ++d)
      X[d] = x00;
  }

  // Constructors for fixed dimension
  inline TinyVector(const T& x00, const T& x01)
  {
    X[0] = x00;
    X[1] = x01;
  }
  inline TinyVector(const T& x00, const T& x01, const T& x02)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
  }
  inline TinyVector(const T& x00, const T& x01, const T& x02, const T& x03)
  {
    X[0] = x00;
    X[1] = x01;
    X[2] = x02;
    X[3] = x03;
  }

  inline TinyVector(const T* restrict base, int offset)
  {
#pragma unroll
    for (int i = 0; i < D; ++i)
      X[i] = base[i * offset];
  }

  inline int size() const { return D; }

  inline TinyVector& operator=(const TinyVector& rhs) = default;
  inline TinyVector& operator=(TinyVector&& rhs)      = default;

  template<class T1>
  inline TinyVector<T, D>& operator=(const TinyVector<T1, D>& rhs)
  {
    for (size_t d = 0; d < D; ++d)
      X[d] = rhs[d];
    return *this;
  }

  inline TinyVector<T, D>& operator=(const T& rhs)
  {
    for (size_t d = 0; d < D; ++d)
      X[d] = rhs;
    return *this;
  }

  // Get and Set Operations
  inline Type_t& operator[](unsigned int i) { return X[i]; }
  inline const Type_t& operator[](unsigned int i) const { return X[i]; }
  inline Type_t* data() { return X; }
  inline const Type_t* data() const { return X; }
  inline Type_t* begin() { return X; }
  inline const Type_t* begin() const { return X; }
  inline Type_t* end() { return X + D; }
  inline const Type_t* end() const { return X + D; }

  TinyVector operator-() const
  {
    TinyVector inverse;
    for (size_t d = 0; d < D; ++d)
      inverse[d] = -X[d];
    return inverse;
  }

  // Elementwise comparison
  bool operator==(const TinyVector<T, D>& that) const
  {
    for (int i = 0; i < D; ++i)
    {
      if ((*this)[i] != that[i])
        return false;
    }
    return true;
  }

  bool operator!=(const TinyVector<T, D>& that) const
  {
    for (int i = 0; i < D; ++i)
    {
      if ((*this)[i] == that[i])
        return false;
    }
    return true;
  }
};

// Adding binary operators using macro defined in OhmmsTinyMeta.h
OHMMS_META_ACCUM_OPERATORS(TinyVector, operator+=, OpAddAssign)
OHMMS_META_ACCUM_OPERATORS(TinyVector, operator-=, OpSubtractAssign)
OHMMS_META_ACCUM_OPERATORS(TinyVector, operator*=, OpMultiplyAssign)
OHMMS_META_ACCUM_OPERATORS(TinyVector, operator/=, OpDivideAssign)

OHMMS_META_BINARY_OPERATORS(TinyVector, operator+, OpAdd)
OHMMS_META_BINARY_OPERATORS(TinyVector, operator-, OpSubtract)
OHMMS_META_BINARY_OPERATORS(TinyVector, operator*, OpMultiply)
OHMMS_META_BINARY_OPERATORS(TinyVector, operator/, OpDivide)

//----------------------------------------------------------------------
// dot product
//----------------------------------------------------------------------
template<class T1, class T2, unsigned D>
inline typename BinaryReturn<T1, T2, OpMultiply>::Type_t dot(const TinyVector<T1, D>& lhs, const TinyVector<T2, D>& rhs)
{
  return OTDot<TinyVector<T1, D>, TinyVector<T2, D>>::apply(lhs, rhs);
}

//----------------------------------------------------------------------
// cross product
//----------------------------------------------------------------------

template<class T1, class T2, unsigned D>
inline TinyVector<typename BinaryReturn<T1, T2, OpMultiply>::Type_t, D> cross(const TinyVector<T1, D>& lhs,
                                                                              const TinyVector<T2, D>& rhs)
{
  return OTCross<TinyVector<T1, D>, TinyVector<T2, D>>::apply(lhs, rhs);
}

//----------------------------------------------------------------------
// cross product
//----------------------------------------------------------------------

template<class T1, class T2, unsigned D>
inline Tensor<typename BinaryReturn<T1, T2, OpMultiply>::Type_t, D> outerProduct(const TinyVector<T1, D>& lhs,
                                                                                 const TinyVector<T2, D>& rhs)
{
  return OuterProduct<TinyVector<T1, D>, TinyVector<T2, D>>::apply(lhs, rhs);
}

//----------------------------------------------------------------------
// I/O
template<class T>
struct printTinyVector
{};

// specialized for Vector<TinyVector<T,D> >
template<class T, unsigned D>
struct printTinyVector<TinyVector<T, D>>
{
  inline static void print(std::ostream& os, const TinyVector<T, D>& r)
  {
    for (int d = 0; d < D; d++)
      os << std::setw(18) << std::setprecision(10) << r[d];
  }
};

template<class T, unsigned D>
std::ostream& operator<<(std::ostream& out, const TinyVector<T, D>& rhs)
{
  printTinyVector<TinyVector<T, D>>::print(out, rhs);
  return out;
}

template<class T, unsigned D>
std::istream& operator>>(std::istream& is, TinyVector<T, D>& rhs)
{
  //printTinyVector<TinyVector<T,D> >::print(out,rhs);
  for (int i = 0; i < D; i++)
    is >> rhs[i];
  return is;
}
} // namespace qmcplusplus

#endif // OHMMS_TINYVECTOR_H
