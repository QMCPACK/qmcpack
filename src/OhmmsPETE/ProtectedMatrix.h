//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_PROTECTEDMATRIX_H
#define QMCPLUSPLUS_PROTECTEDMATRIX_H

#include "OhmmsPETE/OhmmsMatrix.h"

#include <type_traits>

namespace qmcplusplus
{
/** Matrix that cannot be resized after construction
 *
 *  Used for debugging dynamic matrix versus static registered PooledMemory issue
 */
template<class T, typename Alloc = std::allocator<T>>
class ProtectedMatrix : public Matrix<T, Alloc>
{
public:
  using Container = Vector<T, Alloc>;
  using SizeType = typename Container::size_type;
  using Base = Matrix<T, Alloc>;
  ProtectedMatrix() = delete;
  ProtectedMatrix(SizeType n, SizeType m)
  {
    static_assert(std::is_same<T, typename Alloc::value_type>::value, "Matrix and Alloc data types must agree!");
    Base::D1      = n;
    Base::D2      = m;
    Base::TotSize = n * m;
    Base::X.resize(n * m);
  }
  //ProtectedMatrix(const Matrix<T, Alloc>&) { throw std::runtime_error("PProtected Matrix may not be resized"); }
  void resize(SizeType n, SizeType m) { throw std::runtime_error("resize called, Protected Matrix may not be resized"); }
  void add(SizeType n) { throw std::runtime_error("add called, Protected Matrix may not be resized"); }
  void copy(const Matrix<T, Alloc>& rhs) { throw std::runtime_error("copy from Matrix called, Protected Matrix may not be resized"); }
  void copy(const ProtectedMatrix& rhs) { throw std::runtime_error("copy from Protected Matrix called, Protected Matrix may not be resized"); }
  ProtectedMatrix& operator=(const ProtectedMatrix& rhs) { throw std::runtime_error("operator= called, Protected Matrix may not be resized"); }

  template<class RHS, typename Allocator = Alloc, typename = IsHostSafe<Allocator>>
  ProtectedMatrix& operator=(const RHS& rhs) { throw std::runtime_error("templated operator= called, Protected Matrix may not be resized"); }
};

template<class T, typename Alloc>
struct CreateLeaf<ProtectedMatrix<T, Alloc>>
{
  typedef Reference<ProtectedMatrix<T, Alloc>> Leaf_t;
  inline static Leaf_t make(const ProtectedMatrix<T, Alloc>& a) { return Leaf_t(a); }
};

template<class T, typename Alloc>
struct LeafFunctor<ProtectedMatrix<T, Alloc>, EvalLeaf2>
{
  typedef T Type_t;
  inline static Type_t apply(const ProtectedMatrix<T, Alloc>& mat, const EvalLeaf2& f) { return mat(f.val1(), f.val2()); }
};

template<class T, typename Alloc>
struct LeafFunctor<ProtectedMatrix<T, Alloc>, SizeLeaf2>
{
  typedef bool Type_t;
  inline static bool apply(const ProtectedMatrix<T, Alloc>& v, const SizeLeaf2& s) { return s(v.rows(), v.cols()); }
};

}


#endif /* QMCPLUSPLUS_PROTECTEDMATRIX_H */
